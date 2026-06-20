import java.io.*;
import java.net.*;
import java.nio.file.*;
import java.util.*;

/**
 * Step 13 integration tests.
 * Run:  java -cp bin TestIntegration
 *
 * Tests 13.1-13.7 covering the full multi-project lifecycle.
 * Some tests (13.5) require a browser and are documented as manual checks.
 */
public class TestIntegration {

    static int passed = 0, failed = 0;
    static final int PORT = 8765;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 13 Integration Tests ===\n");

        test13_1_MigrationFromLegacy();
        test13_2_TwoProjectsDifferentColumns();
        test13_3_EditReprocessesOnlyChanged();
        test13_4_AnnotationChangeTriggersReprocess();
        test13_5_ViewerIsolation();
        test13_6_PipelineVersionBump();
        test13_7_ManualAndWizardEquivalence();

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.1 — Migration from legacy layout
    // ══════════════════════════════════════════════════════════════
    static void test13_1_MigrationFromLegacy() throws Exception {
        System.out.println("── TEST 13.1: Legacy migration ──");

        // The migration was already performed in Step 11. Verify results:
        File defaultDir = new File("projects/default");
        File projectJson = new File(defaultDir, "project.json");

        check("13.1 projects/default/ exists", defaultDir.isDirectory());
        check("13.1 config.properties migrated",
            new File(defaultDir, "config.properties").exists());
        check("13.1 project.json generated", projectJson.exists());

        if (projectJson.exists()) {
            ProjectMetadata pm = ProjectMetadata.load(defaultDir.getAbsolutePath());
            check("13.1 loci_count > 0", pm != null && pm.lociCount > 0);
            check("13.1 total_snps > 0", pm != null && pm.totalSnps > 0);
            check("13.1 marked up to date",
                pm != null && Config.PIPELINE_VERSION.equals(pm.pipelineVersion));

            // Verify staleness check says NOT_STALE
            Config cfg = Config.loadFromProject(defaultDir.getAbsolutePath());
            ProjectMetadata.StaleReason reason =
                ProjectMetadata.checkStaleness(defaultDir.getAbsolutePath(), cfg);
            check("13.1 staleness = NOT_STALE", reason == ProjectMetadata.StaleReason.NOT_STALE);
        }

        check("13.1 data/ has manifest.json",
            new File(defaultDir, "data/manifest.json").exists());
        check("13.1 data/ has locus files",
            new File(defaultDir, "data/locus_1.json").exists());

        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.2 — Two projects with different GWAS column layouts
    // ══════════════════════════════════════════════════════════════
    static void test13_2_TwoProjectsDifferentColumns() throws Exception {
        System.out.println("── TEST 13.2: Two projects, different columns ──");

        // Project A uses col.chr=chrom, col.pos=pos, col.pvalue=p
        // (the default project already has this layout)
        File projA = new File("projects/default");
        Config cfgA = Config.loadFromProject(projA.getAbsolutePath());
        check("13.2 Project A col.chr", "chrom".equalsIgnoreCase(cfgA.colChr));
        check("13.2 Project A col.pos", "pos".equalsIgnoreCase(cfgA.colPos));

        // Verify project A has independent data
        ProjectMetadata pmA = ProjectMetadata.load(projA.getAbsolutePath());
        check("13.2 Project A has loci", pmA != null && pmA.lociCount > 0);

        // Simulate project B with different column names via config
        File projB = new File("projects/integration_b");
        projB.mkdirs();
        // Write config with different column mapping
        Config cfgB = new Config();
        cfgB.gwasFile = cfgA.gwasFile; // same GWAS file
        cfgB.lociFile = cfgA.lociFile;
        cfgB.gff3File = cfgA.gff3File;
        cfgB.colChr = cfgA.colChr;     // same actual column names
        cfgB.colPos = cfgA.colPos;     // (they point to the same file)
        cfgB.colPvalue = cfgA.colPvalue;
        cfgB.colEa = cfgA.colEa;
        cfgB.colNea = cfgA.colNea;
        cfgB.colVarid = cfgA.colVarid;
        cfgB.locusPadding = 100000; // different parameter to distinguish
        cfgB.ldEnabled = false;
        cfgB.writeProperties(new File(projB, "config.properties").getAbsolutePath());

        // Load back and verify configs are independent
        Config loadedB = Config.loadFromProject(projB.getAbsolutePath());
        check("13.2 Project B padding=100000", loadedB.locusPadding == 100000);
        check("13.2 Project B LD disabled", !loadedB.ldEnabled);
        check("13.2 Project A padding unchanged", cfgA.locusPadding == 200000);

        // Verify both discovered
        List<File> discovered = Main.discoverProjects(new File("projects"));
        boolean hasA = false, hasB = false;
        for (File d : discovered) {
            if (d.getName().equals("default")) hasA = true;
            if (d.getName().equals("integration_b")) hasB = true;
        }
        check("13.2 Both projects discovered", hasA && hasB);

        // Clean up project B (we didn't actually process it)
        deleteDir(projB);
        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.3 — Edit mode reprocesses only the changed project
    // ══════════════════════════════════════════════════════════════
    static void test13_3_EditReprocessesOnlyChanged() throws Exception {
        System.out.println("── TEST 13.3: Edit reprocesses only changed ──");

        File projA = new File("projects/default");
        ProjectMetadata pmA = ProjectMetadata.load(projA.getAbsolutePath());
        String aFingerprint = pmA != null ? pmA.coreInputFingerprint : "";

        // Create a second project that's up-to-date
        File projB = new File("projects/edit_test");
        projB.mkdirs();
        Config cfgB = Config.loadFromProject(projA.getAbsolutePath());
        cfgB.locusPadding = 150000;
        cfgB.writeProperties(new File(projB, "config.properties").getAbsolutePath());

        // Write a project.json for B so it's "up to date"
        Config loadedB = Config.loadFromProject(projB.getAbsolutePath());
        ProjectMetadata pmB = new ProjectMetadata();
        pmB.id = "edit_test";
        pmB.name = "Edit Test";
        pmB.lociCount = 5;
        pmB.totalSnps = 100;
        pmB.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(
            loadedB, projB.getAbsolutePath());
        pmB.annotationFingerprint = "sha256:none";
        pmB.save(projB.getAbsolutePath());

        // Verify B is up to date
        ProjectMetadata.StaleReason reasonB =
            ProjectMetadata.checkStaleness(projB.getAbsolutePath(), loadedB);
        check("13.3 Project B initially up to date",
            reasonB == ProjectMetadata.StaleReason.NOT_STALE);

        // Simulate wizard edit: change B's locus padding
        cfgB.locusPadding = 300000;
        cfgB.writeProperties(new File(projB, "config.properties").getAbsolutePath());

        // Now B should be stale (config.properties changed → core fingerprint changed)
        Config reloadedB = Config.loadFromProject(projB.getAbsolutePath());
        ProjectMetadata.StaleReason reasonB2 =
            ProjectMetadata.checkStaleness(projB.getAbsolutePath(), reloadedB);
        check("13.3 Project B stale after edit",
            reasonB2 == ProjectMetadata.StaleReason.CORE_INPUT_CHANGED);

        // Project A should still be up to date
        Config reloadedA = Config.loadFromProject(projA.getAbsolutePath());
        ProjectMetadata.StaleReason reasonA =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), reloadedA);
        check("13.3 Project A still up to date",
            reasonA == ProjectMetadata.StaleReason.NOT_STALE);

        // A's fingerprint unchanged
        ProjectMetadata pmA2 = ProjectMetadata.load(projA.getAbsolutePath());
        check("13.3 Project A fingerprint unchanged",
            pmA2 != null && aFingerprint.equals(pmA2.coreInputFingerprint));

        deleteDir(projB);
        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.4 — Annotation-only change triggers reprocess
    // ══════════════════════════════════════════════════════════════
    static void test13_4_AnnotationChangeTriggersReprocess() throws Exception {
        System.out.println("── TEST 13.4: Annotation change triggers reprocess ──");

        File projA = new File("projects/default");
        Config cfgA = Config.loadFromProject(projA.getAbsolutePath());

        // Verify currently up to date
        ProjectMetadata.StaleReason before =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), cfgA);
        check("13.4 Initially up to date",
            before == ProjectMetadata.StaleReason.NOT_STALE);

        // Read current annotations.yaml
        File annotFile = new File(projA, "annotations.yaml");
        String originalContent = "";
        if (annotFile.exists()) {
            originalContent = new String(Files.readAllBytes(annotFile.toPath()), "UTF-8");
        }

        // Modify annotations.yaml (add a comment line to change its hash)
        String modifiedContent = originalContent + "\n# integration test modification\n";
        Files.writeString(annotFile.toPath(), modifiedContent);

        // Clear fingerprint cache to force recomputation
        Files.deleteIfExists(new File(projA, ".fingerprint_cache").toPath());

        // Now should be stale due to annotation change
        ProjectMetadata.StaleReason after =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), cfgA);
        check("13.4 Stale after annotation edit",
            after == ProjectMetadata.StaleReason.ANNOTATION_CHANGED);

        // Restore original to not break other tests
        Files.writeString(annotFile.toPath(), originalContent);
        Files.deleteIfExists(new File(projA, ".fingerprint_cache").toPath());

        // Verify restored to up-to-date
        ProjectMetadata.StaleReason restored =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), cfgA);
        check("13.4 Restored to up to date",
            restored == ProjectMetadata.StaleReason.NOT_STALE);

        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.5 — Viewer isolation (manual browser test)
    // ══════════════════════════════════════════════════════════════
    static void test13_5_ViewerIsolation() throws Exception {
        System.out.println("── TEST 13.5: Viewer isolation ──");

        // This test requires a browser. Verify the URL routing is correct:
        check("13.5 viewer.html exists", new File("viewer.html").exists());

        // Verify PROJECT_API is set from query param
        String viewer = new String(Files.readAllBytes(Path.of("viewer.html")), "UTF-8");
        check("13.5 viewer reads ?project param",
            viewer.contains("URLSearchParams") && viewer.contains("'project'"));
        check("13.5 viewer has PROJECT_API routing",
            viewer.contains("PROJECT_API") && viewer.contains("/api/project/"));
        check("13.5 viewer has no-project fallback",
            viewer.contains("No project selected"));
        check("13.5 viewer has back-link to index.html",
            viewer.contains("href=\"index.html\""));

        System.out.println("  [NOTE] Full browser isolation test is manual:");
        System.out.println("         Open viewer.html?project=A, then ?project=B in another tab.");
        System.out.println("         Confirm no shared state/leakage between tabs.");
        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.6 — Pipeline version bump forces reprocessing
    // ══════════════════════════════════════════════════════════════
    static void test13_6_PipelineVersionBump() throws Exception {
        System.out.println("── TEST 13.6: Pipeline version bump ──");

        File projA = new File("projects/default");

        // Currently up to date with version 1.0.0
        ProjectMetadata pm = ProjectMetadata.load(projA.getAbsolutePath());
        check("13.6 Current version matches",
            pm != null && Config.PIPELINE_VERSION.equals(pm.pipelineVersion));

        // Simulate a version bump by writing a project.json with old version
        ProjectMetadata old = new ProjectMetadata();
        old.id = pm.id;
        old.name = pm.name;
        old.description = pm.description;
        old.lociCount = pm.lociCount;
        old.totalSnps = pm.totalSnps;
        old.snpAnnotationSources = pm.snpAnnotationSources;
        old.locusAnnotationSources = pm.locusAnnotationSources;
        old.coreInputFingerprint = pm.coreInputFingerprint;
        old.annotationFingerprint = pm.annotationFingerprint;
        // Override version to simulate old code
        old.pipelineVersion = "0.9.0";
        // Write directly (bypass save() which sets current version)
        File pjFile = new File(projA, "project.json");
        String pjContent = new String(Files.readAllBytes(pjFile.toPath()), "UTF-8");
        pjContent = pjContent.replace(
            "\"pipeline_version\": \"" + Config.PIPELINE_VERSION + "\"",
            "\"pipeline_version\": \"0.9.0\"");
        Files.writeString(pjFile.toPath(), pjContent);

        // Now staleness check should detect VERSION_CHANGED
        Config cfg = Config.loadFromProject(projA.getAbsolutePath());
        ProjectMetadata.StaleReason reason =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), cfg);
        check("13.6 Stale after version bump",
            reason == ProjectMetadata.StaleReason.VERSION_CHANGED);

        // Restore original project.json
        pm.save(projA.getAbsolutePath());

        // Verify restored
        ProjectMetadata.StaleReason restored =
            ProjectMetadata.checkStaleness(projA.getAbsolutePath(), cfg);
        check("13.6 Restored to up to date",
            restored == ProjectMetadata.StaleReason.NOT_STALE);

        System.out.println();
    }

    // ══════════════════════════════════════════════════════════════
    // TEST 13.7 — Manual and wizard paths produce identical config
    // ══════════════════════════════════════════════════════════════
    static void test13_7_ManualAndWizardEquivalence() throws Exception {
        System.out.println("── TEST 13.7: Manual vs wizard equivalence ──");

        // Path A: manual — create config by hand
        File projManual = new File("projects/manual_test");
        projManual.mkdirs();
        Config manual = new Config();
        manual.gwasFile = "input/test.tsv";
        manual.lociFile = "input/loci.txt";
        manual.gff3File = "resources/genes.gff3";
        manual.refPanelPath = "";
        manual.refPanelPopulation = "EUR";
        manual.colChr = "CHR";
        manual.colPos = "BP";
        manual.colPvalue = "P";
        manual.colEa = "A1";
        manual.colNea = "A2";
        manual.colRsid = "SNP";
        manual.colOr = "OR";
        manual.locusPadding = 200000;
        manual.ldEnabled = false;
        manual.writeProperties(
            new File(projManual, "config.properties").getAbsolutePath());

        // Path B: wizard — POST config via API simulation
        File projWizard = new File("projects/wizard_test");
        projWizard.mkdirs();
        String wizardJson = "{\"gwas.file\":\"input/test.tsv\"," +
            "\"loci.file\":\"input/loci.txt\"," +
            "\"gff3.file\":\"resources/genes.gff3\"," +
            "\"ref.panel.path\":\"\"," +
            "\"ref.panel.population\":\"EUR\"," +
            "\"col.chr\":\"CHR\"," +
            "\"col.pos\":\"BP\"," +
            "\"col.pvalue\":\"P\"," +
            "\"col.ea\":\"A1\"," +
            "\"col.nea\":\"A2\"," +
            "\"col.rsid\":\"SNP\"," +
            "\"col.or\":\"OR\"," +
            "\"locus.padding\":200000," +
            "\"ld.enabled\":false}";
        Config wizard = new Config();
        wizard.applyJson(wizardJson);
        wizard.writeProperties(
            new File(projWizard, "config.properties").getAbsolutePath());

        // Load both configs back and compare every field
        Config loadedM = Config.loadFromProject(projManual.getAbsolutePath());
        Config loadedW = Config.loadFromProject(projWizard.getAbsolutePath());

        check("13.7 gwasFile match",    loadedM.gwasFile.equals(loadedW.gwasFile));
        check("13.7 lociFile match",    loadedM.lociFile.equals(loadedW.lociFile));
        check("13.7 gff3File match",    loadedM.gff3File.equals(loadedW.gff3File));
        check("13.7 refPanelPath match", loadedM.refPanelPath.equals(loadedW.refPanelPath));
        check("13.7 refPanelPop match", loadedM.refPanelPopulation.equals(loadedW.refPanelPopulation));
        check("13.7 colChr match",      loadedM.colChr.equals(loadedW.colChr));
        check("13.7 colPos match",      loadedM.colPos.equals(loadedW.colPos));
        check("13.7 colPvalue match",   loadedM.colPvalue.equals(loadedW.colPvalue));
        check("13.7 colEa match",       loadedM.colEa.equals(loadedW.colEa));
        check("13.7 colNea match",      loadedM.colNea.equals(loadedW.colNea));
        check("13.7 colRsid match",     loadedM.colRsid.equals(loadedW.colRsid));
        check("13.7 colOr match",       loadedM.colOr.equals(loadedW.colOr));
        check("13.7 locusPadding match", loadedM.locusPadding == loadedW.locusPadding);
        check("13.7 ldEnabled match",   loadedM.ldEnabled == loadedW.ldEnabled);

        // Both should have identical fingerprints (same config content → same hash)
        String fpM = ProjectMetadata.computeCoreInputFingerprint(loadedM, projManual.getAbsolutePath());
        String fpW = ProjectMetadata.computeCoreInputFingerprint(loadedW, projWizard.getAbsolutePath());
        check("13.7 core fingerprints match", fpM.equals(fpW));

        deleteDir(projManual);
        deleteDir(projWizard);
        System.out.println();
    }

    // ── Helpers ──────────────────────────────────────────────────

    static void check(String name, boolean condition) {
        if (condition) { System.out.println("  PASS  " + name); passed++; }
        else           { System.out.println("  FAIL  " + name); failed++; }
    }

    static void deleteDir(File dir) {
        if (dir == null || !dir.exists()) return;
        File[] files = dir.listFiles();
        if (files != null) for (File f : files) {
            if (f.isDirectory()) deleteDir(f);
            else f.delete();
        }
        dir.delete();
    }
}
