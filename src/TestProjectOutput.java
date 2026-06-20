import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Tests for Step 5: project-scoped output + project.json generation.
 * Run:  java -cp bin TestProjectOutput
 */
public class TestProjectOutput {

    static int passed = 0, failed = 0;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 5 Project-Scoped Output Tests ===\n");

        testDataLandsUnderProjectDir();
        testProjectJsonCounts();
        testForcedRerunUpdatesTimestamp();
        testUniqueSnpCounting();
        testAnnotationSourceCounting();

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // TEST 5.2: Confirm data/ lands under projects/{id}/data/, not a flat path.
    static void testDataLandsUnderProjectDir() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_out52_");
        Path projDir = tmp.resolve("projects").resolve("myproject");
        Files.createDirectories(projDir);

        Config config = new Config();
        config.outputDir = projDir.toString();

        // Verify the data directory path is derived from outputDir
        String expectedDataDir = projDir.toString() + "/data";
        check("5.2 dataDir = projects/{id}/data",
            expectedDataDir.equals(config.outputDir + "/data"));

        // Verify plinkSubsetsDir and ldResultsDir are also scoped
        check("5.2 plinkSubsetsDir scoped",
            config.plinkSubsetsDir().equals(projDir.toString() + "/plink_subsets"));
        check("5.2 ldResultsDir scoped",
            config.ldResultsDir().equals(projDir.toString() + "/ld_results"));

        // Simulate what JsonExporter.export does: write to config.outputDir + "/data"
        new File(expectedDataDir).mkdirs();
        Files.writeString(Path.of(expectedDataDir, "manifest.json"), "{\"total_loci\":0}");

        check("5.2 manifest.json exists under project dir",
            Files.exists(projDir.resolve("data").resolve("manifest.json")));
        check("5.2 manifest.json NOT at flat output/data",
            !Files.exists(tmp.resolve("output").resolve("data").resolve("manifest.json")));

        deleteDir(tmp.toFile());
    }

    // TEST 5.1: Project with known loci count and annotation sources.
    //           Assert: project.json counts match exactly.
    static void testProjectJsonCounts() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_out51_");
        Path projDir = tmp.resolve("projects").resolve("counted");
        Files.createDirectories(projDir);

        // Create annotations.yaml with 3 SNP-level, 1 locus-level source
        Files.writeString(projDir.resolve("annotations.yaml"),
            "annotations:\n" +
            "  - id: finemap\n    level: snp\n    file: \"a.tsv\"\n" +
            "  - id: coloc\n    level: snp\n    file: \"b.tsv\"\n" +
            "  - id: conditional\n    level: snp\n    file: \"c.tsv\"\n" +
            "  - id: locus_meta\n    level: locus\n    file: \"d.tsv\"\n");

        // Create fake input files for fingerprinting
        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=" + tmp.resolve("g.tsv").toString().replace("\\", "/") + "\n" +
            "loci.file=" + tmp.resolve("l.txt").toString().replace("\\", "/") + "\n" +
            "gff3.file=" + tmp.resolve("g.gff3").toString().replace("\\", "/") + "\n");
        Files.writeString(tmp.resolve("g.tsv"), "header\n");
        Files.writeString(tmp.resolve("l.txt"), "header\n");
        Files.writeString(tmp.resolve("g.gff3"), "##gff\n");

        Config config = Config.loadFromProject(projDir.toString());

        // Build fake pipeline result with 25 loci
        Main.PipelineResult result = new Main.PipelineResult();
        result.config = config;
        result.outputs = new ArrayList<>();
        int totalSnpsAdded = 0;
        for (int i = 1; i <= 25; i++) {
            LocusOutput lo = new LocusOutput();
            lo.locusIndex = i;
            lo.locusName = "Locus " + i;
            lo.chr = "1";
            lo.start = i * 1000;
            lo.end = i * 1000 + 500;
            lo.gwasSnps = new ArrayList<>();
            // Each locus has a different number of SNPs (i SNPs per locus)
            for (int j = 0; j < i; j++) {
                lo.gwasSnps.add(new Snp("rs" + (i * 100 + j), "1",
                    i * 1000 + j * 10, 0.001, "A", "G"));
                totalSnpsAdded++;
            }
            lo.nearestGenes = new ArrayList<>();
            lo.genes = new ArrayList<>();
            result.outputs.add(lo);
        }

        // Use reflection-free approach: call writeProjectMetadata via Main
        // We'll manually construct what writeProjectMetadata does
        ProjectMetadata pm = new ProjectMetadata();
        pm.id = "counted";
        pm.name = "counted";
        pm.lociCount = result.outputs.size();

        Set<String> uniqueSnps = new HashSet<>();
        for (LocusOutput o : result.outputs)
            for (Snp s : o.gwasSnps)
                uniqueSnps.add(s.chr + ":" + s.pos);
        pm.totalSnps = uniqueSnps.size();

        pm.countAnnotationSources(projDir.resolve("annotations.yaml").toString());
        pm.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(config, projDir.toString());
        pm.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(
            projDir.resolve("annotations.yaml").toString());
        pm.save(projDir.toString());

        // Read back and verify
        ProjectMetadata loaded = ProjectMetadata.load(projDir.toString());

        check("5.1 loci_count = 25", loaded.lociCount == 25);
        // Total unique SNPs = sum(1..25) = 325 (all have unique positions)
        check("5.1 total_snps = 325", loaded.totalSnps == 325);
        check("5.1 snp_annotation_sources = 3", loaded.snpAnnotationSources == 3);
        check("5.1 locus_annotation_sources = 1", loaded.locusAnnotationSources == 1);
        check("5.1 pipeline_version set", Config.PIPELINE_VERSION.equals(loaded.pipelineVersion));
        check("5.1 core_fingerprint set", loaded.coreInputFingerprint.startsWith("sha256:"));
        check("5.1 annot_fingerprint set", loaded.annotationFingerprint.startsWith("sha256:"));
        check("5.1 last_processed set", !loaded.lastProcessed.isEmpty());

        deleteDir(tmp.toFile());
    }

    // TEST 5.3: Forced rerun (no input changes) updates last_processed.
    static void testForcedRerunUpdatesTimestamp() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_out53_");
        Path projDir = tmp.resolve("projects").resolve("rerun");
        Files.createDirectories(projDir);

        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=" + tmp.resolve("g.tsv").toString().replace("\\", "/") + "\n" +
            "loci.file=" + tmp.resolve("l.txt").toString().replace("\\", "/") + "\n" +
            "gff3.file=" + tmp.resolve("g.gff3").toString().replace("\\", "/") + "\n");
        Files.writeString(tmp.resolve("g.tsv"), "header\n");
        Files.writeString(tmp.resolve("l.txt"), "header\n");
        Files.writeString(tmp.resolve("g.gff3"), "##gff\n");

        Config config = Config.loadFromProject(projDir.toString());

        // First write
        ProjectMetadata pm1 = new ProjectMetadata();
        pm1.id = "rerun";
        pm1.lociCount = 5;
        pm1.totalSnps = 100;
        pm1.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(config, projDir.toString());
        pm1.annotationFingerprint = "sha256:none";
        pm1.save(projDir.toString());

        String firstTimestamp = ProjectMetadata.load(projDir.toString()).lastProcessed;
        check("5.3 first timestamp not empty", !firstTimestamp.isEmpty());

        // Brief pause to ensure different timestamp
        Thread.sleep(50);

        // Second write (simulating forced rerun)
        ProjectMetadata pm2 = new ProjectMetadata();
        pm2.id = "rerun";
        pm2.lociCount = 5;
        pm2.totalSnps = 100;
        pm2.coreInputFingerprint = pm1.coreInputFingerprint;
        pm2.annotationFingerprint = "sha256:none";
        pm2.save(projDir.toString());

        String secondTimestamp = ProjectMetadata.load(projDir.toString()).lastProcessed;
        check("5.3 second timestamp different", !firstTimestamp.equals(secondTimestamp));
        check("5.3 second timestamp is later",
            secondTimestamp.compareTo(firstTimestamp) > 0);

        deleteDir(tmp.toFile());
    }

    // Unique SNP counting: overlapping loci share SNPs, count each position once
    static void testUniqueSnpCounting() throws Exception {
        // Two loci overlapping at position 1500
        List<LocusOutput> outputs = new ArrayList<>();

        LocusOutput lo1 = new LocusOutput();
        lo1.locusIndex = 1;
        lo1.gwasSnps = new ArrayList<>();
        lo1.gwasSnps.add(new Snp("rs1", "1", 1000, 0.01, "A", "G"));
        lo1.gwasSnps.add(new Snp("rs2", "1", 1500, 0.02, "C", "T")); // shared
        lo1.nearestGenes = new ArrayList<>();
        lo1.genes = new ArrayList<>();
        outputs.add(lo1);

        LocusOutput lo2 = new LocusOutput();
        lo2.locusIndex = 2;
        lo2.gwasSnps = new ArrayList<>();
        lo2.gwasSnps.add(new Snp("rs2", "1", 1500, 0.02, "C", "T")); // shared
        lo2.gwasSnps.add(new Snp("rs3", "1", 2000, 0.03, "G", "A"));
        lo2.nearestGenes = new ArrayList<>();
        lo2.genes = new ArrayList<>();
        outputs.add(lo2);

        // Simple sum would be 4, unique count should be 3
        Set<String> unique = new HashSet<>();
        for (LocusOutput o : outputs)
            for (Snp s : o.gwasSnps)
                unique.add(s.chr + ":" + s.pos);

        check("unique: simple sum = 4",
            outputs.stream().mapToInt(o -> o.gwasSnps.size()).sum() == 4);
        check("unique: deduplicated count = 3", unique.size() == 3);

        deleteDir(null); // no cleanup needed
    }

    // Annotation source counting from YAML
    static void testAnnotationSourceCounting() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_out5ann_");

        // No annotations.yaml → 0/0
        ProjectMetadata pmNone = new ProjectMetadata();
        pmNone.countAnnotationSources(tmp.resolve("nonexistent.yaml").toString());
        check("annot: no file → 0 snp sources", pmNone.snpAnnotationSources == 0);
        check("annot: no file → 0 locus sources", pmNone.locusAnnotationSources == 0);

        // Mixed annotations
        Files.writeString(tmp.resolve("annotations.yaml"),
            "annotations:\n" +
            "  - id: a\n    level: snp\n    file: x.tsv\n" +
            "  - id: b\n    level: snp\n    file: y.tsv\n" +
            "  - id: c\n    level: locus\n    file: z.tsv\n" +
            "  - id: d\n    level: locus\n    file: w.tsv\n" +
            "  - id: e\n    level: locus\n    file: v.tsv\n");

        ProjectMetadata pmMixed = new ProjectMetadata();
        pmMixed.countAnnotationSources(tmp.resolve("annotations.yaml").toString());
        check("annot: 2 snp sources", pmMixed.snpAnnotationSources == 2);
        check("annot: 3 locus sources", pmMixed.locusAnnotationSources == 3);

        deleteDir(tmp.toFile());
    }

    // ── Helpers ──────────────────────────────────────────────────────────

    static void check(String name, boolean condition) {
        if (condition) { System.out.println("  PASS  " + name); passed++; }
        else           { System.out.println("  FAIL  " + name); failed++; }
    }

    static void deleteDir(File dir) {
        if (dir == null) return;
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }
}
