import java.io.*;
import java.nio.file.*;

/**
 * Quick smoke test for ProjectMetadata fingerprinting (Step 1 tests).
 * Run:  java -cp bin TestFingerprint
 */
public class TestFingerprint {

    static int passed = 0, failed = 0;

    public static void main(String[] args) throws Exception {
        // Setup: create a temp project directory with minimal files
        Path tmpDir = Files.createTempDirectory("lynx_fp_test_");
        String projectDir = tmpDir.toString();

        // Create minimal config.properties
        Path configPath = tmpDir.resolve("config.properties");
        Files.writeString(configPath,
            "gwas.file=" + tmpDir.resolve("gwas.tsv") + "\n" +
            "loci.file=" + tmpDir.resolve("loci.txt") + "\n" +
            "gff3.file=" + tmpDir.resolve("genes.gff3") + "\n");

        // Create minimal input files
        Files.writeString(tmpDir.resolve("gwas.tsv"), "chrom\tpos\tp\nX\t100\t0.05\n");
        Files.writeString(tmpDir.resolve("loci.txt"), "meta_chr\tmeta_start\tmeta_end\n1\t100\t200\n");
        Files.writeString(tmpDir.resolve("genes.gff3"), "##gff-version 3\n");

        // Create annotations.yaml referencing a TSV
        Files.writeString(tmpDir.resolve("test_annot.tsv"), "SNP\tPIP\nrs1\t0.9\n");
        Files.writeString(tmpDir.resolve("annotations.yaml"),
            "annotations:\n  - id: test\n    level: snp\n    file: \"test_annot.tsv\"\n");

        Config config = new Config();
        config.gwasFile = tmpDir.resolve("gwas.tsv").toString();
        config.lociFile = tmpDir.resolve("loci.txt").toString();
        config.gff3File = tmpDir.resolve("genes.gff3").toString();
        config.refPanelPath = "";

        System.out.println("=== Step 1 Fingerprint Tests ===\n");

        // TEST 1.1: Same files → identical fingerprint
        String fp1 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        String fp2 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        check("1.1 Identical fingerprint on repeat", fp1.equals(fp2));

        // TEST 1.2: Touch mtime only (same content) → fingerprint unchanged
        Thread.sleep(50);
        File gwasFile = tmpDir.resolve("gwas.tsv").toFile();
        gwasFile.setLastModified(System.currentTimeMillis());
        // Clear cache to force re-evaluation
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));
        String fp3 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        check("1.2 Mtime-only touch → same fingerprint", fp1.equals(fp3));

        // TEST 1.3: Change one byte in loci file → fingerprint changes
        Files.writeString(tmpDir.resolve("loci.txt"), "meta_chr\tmeta_start\tmeta_end\n1\t100\t300\n");
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));
        String fp4 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        check("1.3 Loci file content changed → fingerprint changed", !fp1.equals(fp4));

        // Restore for next tests
        Files.writeString(tmpDir.resolve("loci.txt"), "meta_chr\tmeta_start\tmeta_end\n1\t100\t200\n");
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));

        // TEST 1.4: Edit annotations.yaml → annotation fingerprint changes, core doesn't
        String annotPath = tmpDir.resolve("annotations.yaml").toString();
        String annotFp1 = ProjectMetadata.computeAnnotationFingerprint(annotPath);
        String coreFp1  = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);

        Files.writeString(tmpDir.resolve("annotations.yaml"),
            "annotations:\n  - id: test\n    level: snp\n    file: \"test_annot.tsv\"\n" +
            "  - id: extra\n    level: locus\n    file: \"test_annot.tsv\"\n");
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));

        String annotFp2 = ProjectMetadata.computeAnnotationFingerprint(annotPath);
        String coreFp2  = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        check("1.4a Annotation YAML changed → annotation fingerprint changed", !annotFp1.equals(annotFp2));
        check("1.4b Annotation YAML changed → core fingerprint unchanged", coreFp1.equals(coreFp2));

        // TEST 1.5: Edit referenced TSV, leave YAML untouched → annotation fingerprint changes
        Files.writeString(tmpDir.resolve("annotations.yaml"),
            "annotations:\n  - id: test\n    level: snp\n    file: \"test_annot.tsv\"\n");
        String annotFp3 = ProjectMetadata.computeAnnotationFingerprint(annotPath);
        Files.writeString(tmpDir.resolve("test_annot.tsv"), "SNP\tPIP\nrs1\t0.95\n");
        String annotFp4 = ProjectMetadata.computeAnnotationFingerprint(annotPath);
        check("1.5 Referenced TSV changed → annotation fingerprint changed", !annotFp3.equals(annotFp4));

        // TEST 1.6: Pipeline version check
        ProjectMetadata pm = new ProjectMetadata();
        pm.pipelineVersion = "old_version";
        check("1.6 Old pipeline version → isStale", pm.isStale());
        pm.pipelineVersion = Config.PIPELINE_VERSION;
        check("1.6b Current pipeline version → not stale", !pm.isStale());

        // TEST 1.7: Ref panel size+mtime is fast (no SHA-256 of contents)
        // Create fake large ref panel files (just check they don't get hashed)
        Path fakeBed = tmpDir.resolve("panel.bed");
        Files.writeString(fakeBed, "fake bed data");
        Files.writeString(tmpDir.resolve("panel.bim"), "fake bim data");
        Files.writeString(tmpDir.resolve("panel.fam"), "fake fam data");
        config.refPanelPath = tmpDir.resolve("panel").toString();
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));

        long t0 = System.nanoTime();
        String fpRef1 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        long elapsed1 = System.nanoTime() - t0;

        // Touch ref panel mtime → fingerprint should change (size+mtime tracking)
        Thread.sleep(50);
        fakeBed.toFile().setLastModified(System.currentTimeMillis());
        Files.deleteIfExists(tmpDir.resolve(".fingerprint_cache"));
        String fpRef2 = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        check("1.7a Ref panel mtime changed → fingerprint changed", !fpRef1.equals(fpRef2));
        check("1.7b Fingerprint computed in < 1s", elapsed1 < 1_000_000_000L);

        // TEST: project.json round-trip
        ProjectMetadata meta = new ProjectMetadata();
        meta.id = "test_project";
        meta.name = "Test Project";
        meta.description = "A test";
        meta.lociCount = 25;
        meta.totalSnps = 12345;
        meta.snpAnnotationSources = 3;
        meta.locusAnnotationSources = 1;
        meta.coreInputFingerprint = "sha256:abc123";
        meta.annotationFingerprint = "sha256:def456";
        meta.save(projectDir);

        ProjectMetadata loaded = ProjectMetadata.load(projectDir);
        check("Round-trip: id",          "test_project".equals(loaded.id));
        check("Round-trip: name",        "Test Project".equals(loaded.name));
        check("Round-trip: loci_count",  loaded.lociCount == 25);
        check("Round-trip: total_snps",  loaded.totalSnps == 12345);
        check("Round-trip: snp_sources", loaded.snpAnnotationSources == 3);
        check("Round-trip: core_fp",     "sha256:abc123".equals(loaded.coreInputFingerprint));
        check("Round-trip: annot_fp",    "sha256:def456".equals(loaded.annotationFingerprint));
        check("Round-trip: version",     Config.PIPELINE_VERSION.equals(loaded.pipelineVersion));

        // TEST: annotation source counting
        Files.writeString(tmpDir.resolve("annotations.yaml"),
            "annotations:\n" +
            "  - id: a\n    level: snp\n    file: \"x.tsv\"\n" +
            "  - id: b\n    level: snp\n    file: \"y.tsv\"\n" +
            "  - id: c\n    level: locus\n    file: \"z.tsv\"\n");
        ProjectMetadata countTest = new ProjectMetadata();
        countTest.countAnnotationSources(tmpDir.resolve("annotations.yaml").toString());
        check("Annotation count: 2 snp sources", countTest.snpAnnotationSources == 2);
        check("Annotation count: 1 locus source", countTest.locusAnnotationSources == 1);

        // TEST: staleness detection
        ProjectMetadata fresh = new ProjectMetadata();
        fresh.id = "test";
        fresh.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(config, projectDir);
        fresh.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(
            tmpDir.resolve("annotations.yaml").toString());
        fresh.save(projectDir);

        ProjectMetadata.StaleReason reason = ProjectMetadata.checkStaleness(projectDir, config);
        check("Staleness: fresh project is NOT_STALE", reason == ProjectMetadata.StaleReason.NOT_STALE);

        // Delete project.json → NO_METADATA
        Files.delete(tmpDir.resolve("project.json"));
        reason = ProjectMetadata.checkStaleness(projectDir, config);
        check("Staleness: missing project.json → NO_METADATA", reason == ProjectMetadata.StaleReason.NO_METADATA);

        // Cleanup
        deleteDir(tmpDir.toFile());

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    static void check(String name, boolean condition) {
        if (condition) {
            System.out.println("  PASS  " + name);
            passed++;
        } else {
            System.out.println("  FAIL  " + name);
            failed++;
        }
    }

    static void deleteDir(File dir) {
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }
}
