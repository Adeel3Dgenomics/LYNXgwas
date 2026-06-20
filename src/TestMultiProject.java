import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Tests for Step 2: multi-project discovery, staleness, and orchestration.
 * Run:  java -cp bin TestMultiProject
 */
public class TestMultiProject {

    static int passed = 0, failed = 0;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 2 Multi-Project Tests ===\n");

        testDiscovery();
        testStalenessFiltering();
        testForceFlags();
        testConfigLoadFromProject();
        testNoConfigIgnored();

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // TEST 2.5: A folder under projects/ has no config.properties → ignored
    static void testNoConfigIgnored() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_mp_");
        Path projects = tmp.resolve("projects");
        Files.createDirectories(projects);

        // Create a folder with no config → should be ignored
        Files.createDirectories(projects.resolve("empty_folder"));

        // Create a valid project
        Path projA = projects.resolve("alpha");
        Files.createDirectories(projA);
        writeMinimalConfig(projA);

        List<File> discovered = Main.discoverProjects(projects.toFile());
        check("2.5 Folder without config.properties is ignored",
            discovered.size() == 1 && discovered.get(0).getName().equals("alpha"));

        deleteDir(tmp.toFile());
    }

    // TEST 2.1-like: Discovery finds projects with config.properties
    static void testDiscovery() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_disc_");
        Path projects = tmp.resolve("projects");
        Files.createDirectories(projects);

        // Create 3 project dirs
        for (String name : new String[]{"proj_a", "proj_b", "proj_c"}) {
            Path d = projects.resolve(name);
            Files.createDirectories(d);
            writeMinimalConfig(d);
        }

        // Also create the template file (should be ignored)
        Files.writeString(projects.resolve("config.properties.template"), "# template");

        List<File> discovered = Main.discoverProjects(projects.toFile());
        check("Discovery finds 3 projects", discovered.size() == 3);

        List<String> names = new ArrayList<>();
        for (File f : discovered) names.add(f.getName());
        check("Discovery returns sorted names",
            names.equals(Arrays.asList("proj_a", "proj_b", "proj_c")));

        // Non-existent projects/ dir → empty list
        List<File> none = Main.discoverProjects(new File(tmp.toString(), "nonexistent"));
        check("Discovery on missing dir returns empty", none.isEmpty());

        deleteDir(tmp.toFile());
    }

    // Staleness filtering: fresh project is NOT_STALE, missing metadata is stale
    static void testStalenessFiltering() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_stale_");
        Path projDir = tmp.resolve("projects").resolve("test_proj");
        Files.createDirectories(projDir);

        // Create config + input files
        Path gwas = tmp.resolve("gwas.tsv");
        Path loci = tmp.resolve("loci.txt");
        Path gff3 = tmp.resolve("genes.gff3");
        Files.writeString(gwas, "chrom\tpos\tp\n1\t100\t0.05\n");
        Files.writeString(loci, "meta_chr\tmeta_start\tmeta_end\n1\t100\t200\n");
        Files.writeString(gff3, "##gff-version 3\n");

        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=" + gwas.toString().replace("\\", "/") + "\n" +
            "loci.file=" + loci.toString().replace("\\", "/") + "\n" +
            "gff3.file=" + gff3.toString().replace("\\", "/") + "\n");

        Config config = Config.loadFromProject(projDir.toString());

        // No project.json → stale (NO_METADATA)
        ProjectMetadata.StaleReason r1 = ProjectMetadata.checkStaleness(projDir.toString(), config);
        check("No project.json → NO_METADATA", r1 == ProjectMetadata.StaleReason.NO_METADATA);

        // Write project.json with matching fingerprints → NOT_STALE
        ProjectMetadata pm = new ProjectMetadata();
        pm.id = "test_proj";
        pm.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(config, projDir.toString());
        Files.writeString(projDir.resolve("annotations.yaml"), "annotations: []\n");
        pm.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(
            projDir.resolve("annotations.yaml").toString());
        pm.save(projDir.toString());

        ProjectMetadata.StaleReason r2 = ProjectMetadata.checkStaleness(projDir.toString(), config);
        check("Matching fingerprints → NOT_STALE", r2 == ProjectMetadata.StaleReason.NOT_STALE);

        // Modify an input file → CORE_INPUT_CHANGED
        Files.writeString(loci, "meta_chr\tmeta_start\tmeta_end\n1\t100\t999\n");
        Files.deleteIfExists(projDir.resolve(".fingerprint_cache"));
        ProjectMetadata.StaleReason r3 = ProjectMetadata.checkStaleness(projDir.toString(), config);
        check("Input changed → CORE_INPUT_CHANGED", r3 == ProjectMetadata.StaleReason.CORE_INPUT_CHANGED);

        deleteDir(tmp.toFile());
    }

    // --project= and --all flags
    static void testForceFlags() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_force_");
        Path projects = tmp.resolve("projects");

        // Create two projects
        for (String name : new String[]{"alpha", "beta"}) {
            Path d = projects.resolve(name);
            Files.createDirectories(d);
            writeMinimalConfig(d);
        }

        List<File> all = Main.discoverProjects(projects.toFile());
        check("Force: 2 projects discovered", all.size() == 2);

        // Filter by --project=alpha
        List<File> filtered = new ArrayList<>(all);
        filtered.removeIf(d -> !d.getName().equals("alpha"));
        check("Force: --project=alpha filters to 1", filtered.size() == 1);
        check("Force: filtered project is alpha", filtered.get(0).getName().equals("alpha"));

        // --project=nonexistent → empty
        List<File> empty = new ArrayList<>(all);
        empty.removeIf(d -> !d.getName().equals("nonexistent"));
        check("Force: --project=nonexistent → empty", empty.isEmpty());

        deleteDir(tmp.toFile());
    }

    // Config.loadFromProject sets outputDir to projectDir
    static void testConfigLoadFromProject() throws Exception {
        Path tmp = Files.createTempDirectory("lynx_cfgload_");
        Path projDir = tmp.resolve("myproject");
        Files.createDirectories(projDir);

        Path gwas = tmp.resolve("gwas.tsv");
        Path loci = tmp.resolve("loci.txt");
        Path gff3 = tmp.resolve("genes.gff3");
        Files.writeString(gwas, "chrom\tpos\tp\n1\t100\t0.05\n");
        Files.writeString(loci, "meta_chr\tmeta_start\tmeta_end\n1\t100\t200\n");
        Files.writeString(gff3, "##gff-version 3\n");

        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=" + gwas.toString().replace("\\", "/") + "\n" +
            "loci.file=" + loci.toString().replace("\\", "/") + "\n" +
            "gff3.file=" + gff3.toString().replace("\\", "/") + "\n" +
            "locus.padding=500000\n" +
            "ld.enabled=false\n");

        Config config = Config.loadFromProject(projDir.toString());

        check("loadFromProject: outputDir = projectDir",
            config.outputDir.equals(projDir.toString()));
        check("loadFromProject: gwas.file loaded",
            config.gwasFile.equals(gwas.toString().replace("\\", "/")));
        check("loadFromProject: locus.padding parsed", config.locusPadding == 500000);
        check("loadFromProject: ld.enabled parsed", !config.ldEnabled);

        deleteDir(tmp.toFile());
    }

    // ── Helpers ──────────────────────────────────────────────────────────

    static void writeMinimalConfig(Path projDir) throws IOException {
        // Minimal config that points to nonexistent files (for discovery/staleness tests only)
        Files.writeString(projDir.resolve("config.properties"),
            "gwas.file=input/dummy.tsv\n" +
            "loci.file=input/dummy_loci.txt\n" +
            "gff3.file=resources/dummy.gff3\n");
    }

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
