import java.io.*;
import java.net.*;
import java.nio.file.*;
import java.util.*;

/**
 * Tests for Step 6: project-scoped API endpoints.
 * Run:  java -cp bin TestServerEndpoints
 *
 * Requires: projects/test_demo/ with processed data from earlier test run,
 *           OR creates a minimal test project.
 */
public class TestServerEndpoints {

    static int passed = 0, failed = 0;
    static final int PORT = 8765;

    public static void main(String[] args) throws Exception {
        System.out.println("=== Step 6 Server Endpoint Tests ===\n");

        // Setup: create two minimal test projects
        setupTestProjects();

        LocalServer server = null;
        try {
            server = new LocalServer("output");
            server.start();
            System.out.println("Test server started on port " + PORT + "\n");

            testListProjects();
            testProjectManifest();
            testProjectManifestNoLeakage();
            testConfigRoundTrip();
            testPeekFileHeader();
            testProcessEndpoint();

        } finally {
            if (server != null) server.stop();
            cleanupTestProjects();
        }

        System.out.printf("%n=== Results: %d passed, %d failed ===%n", passed, failed);
        if (failed > 0) System.exit(1);
    }

    // TEST 6.1: GET /api/projects with one stale, one current project
    static void testListProjects() throws Exception {
        String resp = httpGet("/api/projects");
        check("6.1 Response is JSON array", resp.startsWith("["));
        check("6.1 Contains proj_alpha", resp.contains("proj_alpha"));
        check("6.1 Contains proj_beta", resp.contains("proj_beta"));
        // proj_alpha has project.json → up_to_date, proj_beta has none → needs_reprocessing
        check("6.1 proj_alpha status up_to_date",
            resp.contains("\"id\":\"proj_alpha\"") &&
            resp.contains("\"status\":\"up_to_date\""));
        check("6.1 proj_beta status needs_reprocessing",
            resp.contains("\"id\":\"proj_beta\"") &&
            resp.contains("\"status\":\"needs_reprocessing\""));
    }

    // TEST 6.2: GET /api/project/A/manifest vs /api/project/B/manifest — no cross-leakage
    static void testProjectManifest() throws Exception {
        String resp = httpGet("/api/project/proj_alpha/manifest");
        check("6.2 proj_alpha manifest returns JSON", resp.startsWith("{"));
        check("6.2 proj_alpha manifest has total_loci", resp.contains("\"total_loci\""));
    }

    static void testProjectManifestNoLeakage() throws Exception {
        String respA = httpGet("/api/project/proj_alpha/manifest");
        // proj_beta has no data, should 404
        int statusB = httpStatus("/api/project/proj_beta/manifest");
        check("6.2 proj_beta manifest returns 404 (no data)", statusB == 404);
        // Re-fetch A to confirm it's still correct after B's 404
        String respA2 = httpGet("/api/project/proj_alpha/manifest");
        check("6.2 No cross-project leakage", respA.equals(respA2));
    }

    // TEST 6.3: POST /api/project/A/config, then GET — round-trip
    static void testConfigRoundTrip() throws Exception {
        String postBody = "{\"gwas.file\":\"input/new_gwas.tsv\"," +
            "\"loci.file\":\"input/new_loci.txt\"," +
            "\"gff3.file\":\"resources/new_genes.gff3\"," +
            "\"ref.panel.path\":\"\"," +
            "\"ref.panel.population\":\"AFR\"," +
            "\"col.chr\":\"CHROM\"," +
            "\"col.pos\":\"POS\"," +
            "\"col.pvalue\":\"PVAL\"," +
            "\"col.ea\":\"EA\"," +
            "\"col.nea\":\"NEA\"," +
            "\"col.rsid\":\"RS\"," +
            "\"col.varid\":\"\"," +
            "\"col.or\":\"OR\"," +
            "\"locus.padding\":300000," +
            "\"ld.enabled\":false," +
            "\"ld.triangle.boundary\":150}";

        String postResp = httpPost("/api/project/proj_beta/config", postBody);
        check("6.3 POST config returns ok", postResp.contains("\"ok\":true"));

        String getResp = httpGet("/api/project/proj_beta/config");
        check("6.3 GET config returns JSON", getResp.startsWith("{"));
        check("6.3 round-trip: gwas.file", getResp.contains("input/new_gwas.tsv"));
        check("6.3 round-trip: population AFR", getResp.contains("AFR"));
        check("6.3 round-trip: col.chr=CHROM", getResp.contains("CHROM"));
        check("6.3 round-trip: col.pvalue=PVAL", getResp.contains("PVAL"));
        check("6.3 round-trip: col.or=OR", getResp.contains("\"col.or\":\"OR\""));

        // Verify the underlying file is standard config.properties format
        String configContent = Files.readString(
            Path.of("projects/proj_beta/config.properties"));
        check("6.3 file: gwas.file key present", configContent.contains("gwas.file=input/new_gwas.tsv"));
        check("6.3 file: ref.panel.population=AFR", configContent.contains("ref.panel.population=AFR"));
        check("6.3 file: col.or=OR", configContent.contains("col.or=OR"));
    }

    // TEST 6.4: GET /api/peek-file-header
    static void testPeekFileHeader() throws Exception {
        // Use the test GWAS file we created
        String gwasPath = new File("projects/proj_alpha/test_gwas.tsv").getAbsolutePath();
        String encodedPath = URLEncoder.encode(gwasPath, "UTF-8");
        String resp = httpGet("/api/peek-file-header?path=" + encodedPath);
        check("6.4 peek returns JSON", resp.startsWith("{"));
        check("6.4 peek has columns array", resp.contains("\"columns\""));
        check("6.4 peek contains CHR", resp.contains("CHR"));
        check("6.4 peek contains BP", resp.contains("BP"));
        check("6.4 peek contains P", resp.contains("\"P\""));
    }

    // TEST 6.5: POST /api/projects/process
    static void testProcessEndpoint() throws Exception {
        // Just test that it accepts the request — actual processing would need real data
        String resp = httpPost("/api/projects/process", "{\"id\":\"proj_nonexistent\"}");
        check("6.5 process nonexistent returns 404",
            resp.contains("not found") || resp.contains("error"));

        // Process with valid but not-runnable project (will start and fail in background)
        String resp2 = httpPost("/api/projects/process", "{\"id\":\"proj_beta\"}");
        check("6.5 process valid project accepted", resp2.contains("processing"));
    }

    // ── Setup / teardown ─────────────────────────────────────────────────

    static void setupTestProjects() throws Exception {
        // proj_alpha: has data and project.json (current)
        Path alpha = Path.of("projects/proj_alpha");
        Files.createDirectories(alpha.resolve("data"));

        // Create a minimal GWAS file for peek-file-header test
        Files.writeString(alpha.resolve("test_gwas.tsv"),
            "CHR\tBP\tP\tA1\tA2\tSNP\n1\t1000\t0.05\tA\tG\trs1\n");

        Files.writeString(alpha.resolve("config.properties"),
            "gwas.file=" + alpha.resolve("test_gwas.tsv").toAbsolutePath().toString().replace("\\","/") + "\n" +
            "loci.file=input/loci.txt\n" +
            "gff3.file=resources/gencode.v37.annotation.gff3\n");

        // Write a manifest
        Files.writeString(alpha.resolve("data/manifest.json"),
            "{\"total_loci\":2,\"loci\":[{\"index\":1,\"chr\":\"1\",\"start\":100,\"end\":200}]}");

        // Write project.json with matching fingerprints so it's "up to date"
        Config cfgA = Config.loadFromProject(alpha.toAbsolutePath().toString());
        ProjectMetadata pmA = new ProjectMetadata();
        pmA.id = "proj_alpha"; pmA.name = "Alpha Project";
        pmA.lociCount = 2; pmA.totalSnps = 50;
        pmA.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(
            cfgA, alpha.toAbsolutePath().toString());
        pmA.annotationFingerprint = "sha256:none";
        pmA.save(alpha.toAbsolutePath().toString());

        // proj_beta: has config but no project.json (stale/never processed)
        Path beta = Path.of("projects/proj_beta");
        Files.createDirectories(beta);
        Files.writeString(beta.resolve("config.properties"),
            "gwas.file=input/dummy.tsv\n" +
            "loci.file=input/loci.txt\n" +
            "gff3.file=resources/gencode.v37.annotation.gff3\n");
    }

    static void cleanupTestProjects() {
        deleteDir(new File("projects/proj_alpha"));
        deleteDir(new File("projects/proj_beta"));
        // Don't delete projects/ itself — other tests may have left content
        File projects = new File("projects");
        if (projects.isDirectory()) {
            String[] remaining = projects.list();
            if (remaining != null && remaining.length == 0) projects.delete();
        }
    }

    // ── HTTP helpers ─────────────────────────────────────────────────────

    static String httpGet(String path) throws Exception {
        URL url = new URL("http://localhost:" + PORT + path);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        conn.setConnectTimeout(5000);
        conn.setReadTimeout(5000);
        try (InputStream is = conn.getInputStream()) {
            return new String(is.readAllBytes(), "UTF-8");
        }
    }

    static int httpStatus(String path) throws Exception {
        URL url = new URL("http://localhost:" + PORT + path);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        conn.setConnectTimeout(5000);
        conn.setReadTimeout(5000);
        return conn.getResponseCode();
    }

    static String httpPost(String path, String body) throws Exception {
        URL url = new URL("http://localhost:" + PORT + path);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("POST");
        conn.setDoOutput(true);
        conn.setRequestProperty("Content-Type", "application/json");
        conn.setConnectTimeout(5000);
        conn.setReadTimeout(5000);
        try (OutputStream os = conn.getOutputStream()) {
            os.write(body.getBytes("UTF-8"));
        }
        InputStream is = conn.getResponseCode() < 400
            ? conn.getInputStream() : conn.getErrorStream();
        try { return new String(is.readAllBytes(), "UTF-8"); }
        finally { is.close(); }
    }

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
