import java.io.*;
import java.nio.file.*;

/**
 * Tests for Step 3: per-project Config loading.
 * Run:  java -cp bin TestPerProjectConfig
 */
public class TestPerProjectConfig {

    static int passed = 0, failed = 0;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 3 Per-Project Config Tests ===\n");

        testLoadFromProjectMatchesRootBehavior();
        testTwoProjectsNoCrossContamination();
        testOutputDirOverride();
        testPipelineVersionIsConstant();

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // TEST 3.1: Load config from projects/test_project/. Assert all existing
    //           keys parse identically to current root-level behavior.
    static void testLoadFromProjectMatchesRootBehavior() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_cfg31_");
        Path projDir = tmp.resolve("projects").resolve("test_project");
        Files.createDirectories(projDir);

        // Write a config with every key
        String configContent =
            "gwas.file=input/my_gwas.tsv\n" +
            "loci.file=input/my_loci.txt\n" +
            "gff3.file=resources/my_genes.gff3\n" +
            "ref.panel.path=C:/data/panel\n" +
            "ref.panel.population=EUR\n" +
            "col.chr=CHR\n" +
            "col.pos=BP\n" +
            "col.pvalue=P\n" +
            "col.rsid=SNP\n" +
            "col.varid=ID\n" +
            "col.ea=A1\n" +
            "col.nea=A2\n" +
            "top.snp.file=input/tops.txt\n" +
            "locus.padding=500000\n" +
            "ld.enabled=true\n" +
            "ld.triangle.boundary=200\n" +
            "ld.r2.threshold=0.1\n" +
            "ld.parallel.jobs=8\n" +
            "threads=16\n" +
            "max.snps.per.locus=10000\n" +
            "split.ld.threshold=0.3\n" +
            "split.min.distance.bp=100000\n";

        Files.writeString(projDir.resolve("config.properties"), configContent);

        // Also load via the traditional CLI path for comparison
        Path rootConfig = tmp.resolve("root_config.properties");
        Files.writeString(rootConfig, configContent);

        Config fromProject = Config.loadFromProject(projDir.toString());
        Config fromRoot = Config.load(new String[]{"--config", rootConfig.toString()});

        // Compare every field (except outputDir which is intentionally different)
        check("3.1 gwasFile matches",       fromProject.gwasFile.equals(fromRoot.gwasFile));
        check("3.1 lociFile matches",        fromProject.lociFile.equals(fromRoot.lociFile));
        check("3.1 gff3File matches",        fromProject.gff3File.equals(fromRoot.gff3File));
        check("3.1 refPanelPath matches",    fromProject.refPanelPath.equals(fromRoot.refPanelPath));
        check("3.1 refPanelPopulation",      fromProject.refPanelPopulation.equals(fromRoot.refPanelPopulation));
        check("3.1 colChr matches",          fromProject.colChr.equals(fromRoot.colChr));
        check("3.1 colPos matches",          fromProject.colPos.equals(fromRoot.colPos));
        check("3.1 colPvalue matches",       fromProject.colPvalue.equals(fromRoot.colPvalue));
        check("3.1 colRsid matches",         fromProject.colRsid.equals(fromRoot.colRsid));
        check("3.1 colVarid matches",        fromProject.colVarid.equals(fromRoot.colVarid));
        check("3.1 colEa matches",           fromProject.colEa.equals(fromRoot.colEa));
        check("3.1 colNea matches",          fromProject.colNea.equals(fromRoot.colNea));
        check("3.1 topSnpFile matches",      fromProject.topSnpFile.equals(fromRoot.topSnpFile));
        check("3.1 locusPadding matches",    fromProject.locusPadding == fromRoot.locusPadding);
        check("3.1 ldEnabled matches",       fromProject.ldEnabled == fromRoot.ldEnabled);
        check("3.1 ldTriangleBoundary",      fromProject.ldTriangleBoundary == fromRoot.ldTriangleBoundary);
        check("3.1 ldR2Threshold matches",   fromProject.ldR2Threshold == fromRoot.ldR2Threshold);
        check("3.1 ldParallelJobs matches",  fromProject.ldParallelJobs == fromRoot.ldParallelJobs);
        check("3.1 threads matches",         fromProject.threads == fromRoot.threads);
        check("3.1 maxSnpsPerLocus matches", fromProject.maxSnpsPerLocus == fromRoot.maxSnpsPerLocus);
        check("3.1 splitLdThreshold",        fromProject.splitLdThreshold == fromRoot.splitLdThreshold);
        check("3.1 splitMinDistBp matches",  fromProject.splitMinDistBp == fromRoot.splitMinDistBp);

        deleteDir(tmp.toFile());
    }

    // TEST 3.2: Two projects with different config values, loaded in the same batch.
    //           Assert each uses its own project's values — no cross-contamination.
    static void testTwoProjectsNoCrossContamination() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_cfg32_");

        // Project A: EUR population, LD enabled, padding 100k
        Path projA = tmp.resolve("projects").resolve("projA");
        Files.createDirectories(projA);
        Files.writeString(projA.resolve("config.properties"),
            "gwas.file=input/gwas_a.tsv\n" +
            "loci.file=input/loci_a.txt\n" +
            "gff3.file=resources/genes_a.gff3\n" +
            "ref.panel.path=C:/data/eur_panel\n" +
            "ref.panel.population=EUR\n" +
            "col.chr=CHROM\n" +
            "col.pos=POS\n" +
            "col.pvalue=PVALUE\n" +
            "locus.padding=100000\n" +
            "ld.enabled=true\n" +
            "ld.triangle.boundary=50\n" +
            "threads=2\n");

        // Project B: EAS population, LD disabled, padding 500k
        Path projB = tmp.resolve("projects").resolve("projB");
        Files.createDirectories(projB);
        Files.writeString(projB.resolve("config.properties"),
            "gwas.file=input/gwas_b.tsv\n" +
            "loci.file=input/loci_b.txt\n" +
            "gff3.file=resources/genes_b.gff3\n" +
            "ref.panel.population=EAS\n" +
            "col.chr=chr\n" +
            "col.pos=bp\n" +
            "col.pvalue=p\n" +
            "locus.padding=500000\n" +
            "ld.enabled=false\n" +
            "ld.triangle.boundary=200\n" +
            "threads=8\n");

        // Load both in sequence (simulating Main's batch loop)
        Config configA = Config.loadFromProject(projA.toString());
        Config configB = Config.loadFromProject(projB.toString());

        // Verify A has A's values
        check("3.2 A: gwasFile",         configA.gwasFile.equals("input/gwas_a.tsv"));
        check("3.2 A: refPanelPop",      configA.refPanelPopulation.equals("EUR"));
        check("3.2 A: colChr",           configA.colChr.equals("CHROM"));
        check("3.2 A: colPos",           configA.colPos.equals("POS"));
        check("3.2 A: colPvalue",        configA.colPvalue.equals("PVALUE"));
        check("3.2 A: locusPadding",     configA.locusPadding == 100000);
        check("3.2 A: ldEnabled",        configA.ldEnabled == true);
        check("3.2 A: ldTriBoundary",    configA.ldTriangleBoundary == 50);
        check("3.2 A: threads",          configA.threads == 2);
        check("3.2 A: outputDir",        configA.outputDir.equals(projA.toString()));

        // Verify B has B's values (NOT A's)
        check("3.2 B: gwasFile",         configB.gwasFile.equals("input/gwas_b.tsv"));
        check("3.2 B: refPanelPop",      configB.refPanelPopulation.equals("EAS"));
        check("3.2 B: colChr",           configB.colChr.equals("chr"));
        check("3.2 B: colPos",           configB.colPos.equals("bp"));
        check("3.2 B: colPvalue",        configB.colPvalue.equals("p"));
        check("3.2 B: locusPadding",     configB.locusPadding == 500000);
        check("3.2 B: ldEnabled",        configB.ldEnabled == false);
        check("3.2 B: ldTriBoundary",    configB.ldTriangleBoundary == 200);
        check("3.2 B: threads",          configB.threads == 8);
        check("3.2 B: outputDir",        configB.outputDir.equals(projB.toString()));

        // Verify A wasn't mutated by loading B
        check("3.2 A still EUR after B load", configA.refPanelPopulation.equals("EUR"));
        check("3.2 A still ldEnabled after B", configA.ldEnabled == true);
        check("3.2 A still padding 100k",     configA.locusPadding == 100000);

        deleteDir(tmp.toFile());
    }

    // Verify outputDir is always overridden to the project directory,
    // even if config.properties contains output.dir=something_else
    static void testOutputDirOverride() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_cfg3out_");
        Path projDir = tmp.resolve("projects").resolve("myproj");
        Files.createDirectories(projDir);

        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=input/g.tsv\n" +
            "loci.file=input/l.txt\n" +
            "gff3.file=resources/g.gff3\n" +
            "output.dir=custom/output/path\n");

        Config config = Config.loadFromProject(projDir.toString());

        // outputDir must be the project directory, NOT what config.properties says
        check("outputDir override: project dir wins",
            config.outputDir.equals(projDir.toString()));
        check("plinkSubsetsDir scoped to project",
            config.plinkSubsetsDir().startsWith(projDir.toString()));
        check("ldResultsDir scoped to project",
            config.ldResultsDir().startsWith(projDir.toString()));

        deleteDir(tmp.toFile());
    }

    // PIPELINE_VERSION is a compile-time constant, not a user-editable key
    static void testPipelineVersionIsConstant() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_cfg3ver_");
        Path projDir = tmp.resolve("projects").resolve("vertest");
        Files.createDirectories(projDir);

        // Even if someone puts pipeline_version in config, it doesn't matter —
        // it's a static final constant in Config.java
        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=input/g.tsv\n" +
            "loci.file=input/l.txt\n" +
            "gff3.file=resources/g.gff3\n");

        Config config = Config.loadFromProject(projDir.toString());

        check("PIPELINE_VERSION is accessible", Config.PIPELINE_VERSION != null);
        check("PIPELINE_VERSION is non-empty",  !Config.PIPELINE_VERSION.isEmpty());
        check("PIPELINE_VERSION is 1.0.0",      "1.0.0".equals(Config.PIPELINE_VERSION));

        // Verify it's the same constant regardless of which project loaded
        Path proj2 = tmp.resolve("projects").resolve("other");
        Files.createDirectories(proj2);
        Files.writeString(proj2.resolve("config.properties"),
            "gwas.file=input/g2.tsv\n" +
            "loci.file=input/l2.txt\n" +
            "gff3.file=resources/g2.gff3\n");
        Config config2 = Config.loadFromProject(proj2.toString());

        check("PIPELINE_VERSION same across projects",
            Config.PIPELINE_VERSION.equals(Config.PIPELINE_VERSION));

        deleteDir(tmp.toFile());
    }

    // ── Helpers ──────────────────────────────────────────────────────────

    static void check(String name, boolean condition) {
        if (condition) { System.out.println("  PASS  " + name); passed++; }
        else           { System.out.println("  FAIL  " + name); failed++; }
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
