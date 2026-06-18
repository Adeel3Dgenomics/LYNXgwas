import com.sun.net.httpserver.*;
import javax.swing.*;
import java.io.*;
import java.net.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Minimal HTTP server on localhost:8765.
 *
 *   GET  /manifest      → manifest.json
 *   GET  /locus/{n}     → locus_N.json
 *   GET  /pick-folder   → opens JFileChooser, returns {"path":"..."}
 *   POST /save-pdf      → body: {filename,folder,data(base64 dataURI)}, writes file
 *   GET  /progress      → {"phase":"...","locus":n,"total":n,"pct":n,"done":false}
 *   GET  /              → static files from outputDir
 */
public class LocalServer {

    public static final int PORT = 8765;

    /** Shared progress state; Main updates this during the pipeline. */
    public static final ProgressTracker progress = new ProgressTracker();

    private final HttpServer http;
    private final String     dataDir;
    private final String     outputDir;
    private volatile String  lastFolder;

    // Shared state for locus updates (set by Main after pipeline completes)
    private Config             config;
    private GffParser          gff;
    private List<Locus>        loci;
    private List<LocusOutput>  outputs;

    public void setPipelineState(Config config, GffParser gff, List<Locus> loci, List<LocusOutput> outputs) {
        this.config  = config;
        this.gff     = gff;
        this.loci    = loci;
        this.outputs = outputs;
    }

    public LocalServer(String outputDir) throws IOException {
        this.outputDir  = outputDir;
        this.dataDir    = outputDir + "/data";
        this.lastFolder = outputDir + "/plots";
        new File(this.lastFolder).mkdirs();

        http = HttpServer.create(new InetSocketAddress("localhost", PORT), 32);
        http.createContext("/manifest",       this::manifest);
        http.createContext("/locus/",         this::locus);
        http.createContext("/pick-folder",    this::pickFolder);
        http.createContext("/save-pdf",       this::savePdf);
        http.createContext("/progress",       this::progressEndpoint);
        http.createContext("/api/update-locus",  this::updateLocus);
        http.createContext("/api/create-locus", this::createLocus);
        http.createContext("/api/split-locus",    this::splitLocus);
        http.createContext("/api/validate-split", this::validateSplit);
        http.createContext("/",               this::staticFiles);
        http.setExecutor(Executors.newFixedThreadPool(4));
    }

    public void start() { http.start(); }
    public void stop()  { http.stop(0); }

    // ── Endpoints ──────────────────────────────────────────────────────────

    private void manifest(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        serveJson(ex, dataDir + "/manifest.json");
    }

    private void locus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String n = ex.getRequestURI().getPath().replaceAll("\\D", "");
        serveJson(ex, dataDir + "/locus_" + n + ".json");
    }

    private void pickFolder(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String chosen = folderPicker(lastFolder);
        if (chosen != null) lastFolder = chosen;
        String json = chosen != null
            ? "{\"path\":\"" + escape(chosen) + "\"}"
            : "{\"path\":null}";
        respond(ex, 200, "application/json", json.getBytes());
    }

    private void savePdf(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        String json     = new String(readAll(ex.getRequestBody()), "UTF-8");
        String filename = extractStr(json, "filename");
        String folder   = extractStr(json, "folder");
        String data     = extractStr(json, "data");
        if (filename == null || folder == null || data == null) {
            respond(ex, 400, "application/json", "{\"error\":\"missing fields\"}".getBytes()); return;
        }
        int comma = data.indexOf(',');
        byte[] pdf = Base64.getDecoder().decode(
            (comma >= 0 ? data.substring(comma + 1) : data).replaceAll("\\s", ""));
        new File(folder).mkdirs();
        Files.write(Paths.get(folder, filename), pdf);
        respond(ex, 200, "application/json", "{\"ok\":true}".getBytes());
    }

    private void progressEndpoint(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        String json = String.format(
            "{\"phase\":\"%s\",\"locus\":%d,\"total\":%d,\"pct\":%d,\"done\":%s}",
            escape(progress.phase), progress.locusIndex,
            progress.totalLoci, progress.pct(),
            progress.done ? "true" : "false");
        respond(ex, 200, "application/json", json.getBytes());
    }

    private void staticFiles(HttpExchange ex) throws IOException {
        cors(ex);
        String p = ex.getRequestURI().getPath();
        if (p.equals("/")) p = "/index.html";
        File f = new File(outputDir + p);
        if (!f.exists() || !f.isFile()) {
            respond(ex, 404, "text/plain", ("Not found: " + p).getBytes()); return;
        }
        byte[] bytes = Files.readAllBytes(f.toPath());
        ex.getResponseHeaders().set("Content-Type", guessMime(p));
        ex.sendResponseHeaders(200, bytes.length);
        try (OutputStream os = ex.getResponseBody()) { os.write(bytes); }
    }

    private void updateLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sIdx   = extractStr(body, "locus_index");
        String sStart = extractStr(body, "new_start");
        String sEnd   = extractStr(body, "new_end");

        if (sIdx == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing locus_index, new_start, or new_end\"}".getBytes()); return;
        }

        int locusIndex = Integer.parseInt(sIdx);
        long newStart  = Long.parseLong(sStart);
        long newEnd    = Long.parseLong(sEnd);

        LocusUpdater.UpdateResult ur = LocusUpdater.update(
            locusIndex, newStart, newEnd, config, gff, loci);

        if (ur.ok) {
            byte[] json = Files.readAllBytes(new File(ur.jsonPath).toPath());
            respond(ex, 200, "application/json", json);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escape(ur.error) + "\"}").getBytes());
        }
    }

    private void createLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null || outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        String body  = new String(readAll(ex.getRequestBody()), "UTF-8");
        String chr   = extractStr(body, "chr");
        String sStart = extractStr(body, "start");
        String sEnd   = extractStr(body, "end");
        String name   = extractStr(body, "locus_name");

        if (chr == null || sStart == null || sEnd == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing chr, start, or end\"}".getBytes()); return;
        }

        long start = Long.parseLong(sStart);
        long end   = Long.parseLong(sEnd);

        LocusUpdater.UpdateResult ur = LocusUpdater.create(
            chr, start, end, name, config, gff, loci, outputs);

        if (ur.ok) {
            // Return the updated manifest so the frontend can reload
            byte[] manifestBytes = Files.readAllBytes(
                new File(dataDir + "/manifest.json").toPath());
            respond(ex, 200, "application/json", manifestBytes);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escape(ur.error) + "\"}").getBytes());
        }
    }

    private void splitLocus(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || gff == null || loci == null || outputs == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        int origIndex = Integer.parseInt(sOrigIdx);

        // Parse regions array manually: [{name,start,end}, ...]
        List<LocusUpdater.SplitRegion> regions = new ArrayList<>();
        int ri = body.indexOf("\"regions\"");
        if (ri < 0) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing regions\"}".getBytes()); return;
        }
        // Find array start
        int arrStart = body.indexOf('[', ri);
        int arrEnd   = body.indexOf(']', arrStart);
        if (arrStart < 0 || arrEnd < 0) {
            respond(ex, 400, "application/json",
                "{\"error\":\"invalid regions format\"}".getBytes()); return;
        }
        String arrStr = body.substring(arrStart, arrEnd + 1);
        // Split on "},{" to get individual objects
        String[] parts = arrStr.replaceAll("^\\[\\{", "").replaceAll("\\}\\]$", "").split("\\},\\{");
        for (String part : parts) {
            String obj = "{" + part + "}";
            String name  = extractStr(obj, "name");
            String start = extractStr(obj, "start");
            String end   = extractStr(obj, "end");
            if (start == null || end == null) continue;
            regions.add(new LocusUpdater.SplitRegion(
                name != null ? name : "", Long.parseLong(start), Long.parseLong(end)));
        }

        if (regions.size() < 2) {
            respond(ex, 400, "application/json",
                "{\"error\":\"need at least 2 regions to split\"}".getBytes()); return;
        }

        LocusUpdater.UpdateResult ur = LocusUpdater.split(
            origIndex, regions, config, gff, loci, outputs);

        if (ur.ok) {
            byte[] manifestBytes = Files.readAllBytes(
                new File(dataDir + "/manifest.json").toPath());
            respond(ex, 200, "application/json", manifestBytes);
        } else {
            respond(ex, 500, "application/json",
                ("{\"error\":\"" + escape(ur.error) + "\"}").getBytes());
        }
    }

    private void validateSplit(HttpExchange ex) throws IOException {
        cors(ex); if (preflight(ex)) return;
        if (!"POST".equalsIgnoreCase(ex.getRequestMethod())) {
            respond(ex, 405, "text/plain", "Method Not Allowed".getBytes()); return;
        }
        if (config == null || loci == null) {
            respond(ex, 503, "application/json",
                "{\"error\":\"Pipeline state not available\"}".getBytes()); return;
        }

        String body = new String(readAll(ex.getRequestBody()), "UTF-8");
        String sOrigIdx = extractStr(body, "locus_index");
        if (sOrigIdx == null) {
            respond(ex, 400, "application/json",
                "{\"error\":\"missing locus_index\"}".getBytes()); return;
        }
        int origIndex = Integer.parseInt(sOrigIdx);

        // Parse regions
        List<LocusUpdater.SplitRegion> regions = new ArrayList<>();
        int ri = body.indexOf("\"regions\"");
        if (ri >= 0) {
            int arrStart = body.indexOf('[', ri);
            int arrEnd   = body.indexOf(']', arrStart);
            if (arrStart >= 0 && arrEnd >= 0) {
                String arrStr = body.substring(arrStart, arrEnd + 1);
                String[] parts = arrStr.replaceAll("^\\[\\{", "").replaceAll("\\}\\]$", "").split("\\},\\{");
                for (String part : parts) {
                    String obj = "{" + part + "}";
                    String name  = extractStr(obj, "name");
                    String start = extractStr(obj, "start");
                    String end   = extractStr(obj, "end");
                    if (start == null || end == null) continue;
                    regions.add(new LocusUpdater.SplitRegion(
                        name != null ? name : "", Long.parseLong(start), Long.parseLong(end)));
                }
            }
        }

        LocusUpdater.ValidationResult vr = LocusUpdater.validateSplit(
            origIndex, regions, config, loci);

        // Build JSON response
        StringBuilder json = new StringBuilder();
        json.append("{\"valid\":").append(vr.valid);
        json.append(",\"split_ld_threshold\":").append(config.splitLdThreshold);
        json.append(",\"split_min_distance_bp\":").append(config.splitMinDistBp);
        json.append(",\"warnings\":[");
        for (int i = 0; i < vr.warnings.size(); i++) {
            if (i > 0) json.append(',');
            json.append('"').append(escape(vr.warnings.get(i))).append('"');
        }
        json.append("],\"ld_violations\":[");
        for (int i = 0; i < vr.ldViolations.size(); i++) {
            if (i > 0) json.append(',');
            LocusUpdater.LdViolation v = vr.ldViolations.get(i);
            json.append("{\"snp_a\":\"").append(escape(v.snpA))
                .append("\",\"snp_b\":\"").append(escape(v.snpB))
                .append("\",\"pos_a\":").append(v.posA)
                .append(",\"pos_b\":").append(v.posB)
                .append(",\"region_a\":").append(v.regionA)
                .append(",\"region_b\":").append(v.regionB)
                .append(",\"r2\":").append(String.format("%.4f", v.r2))
                .append('}');
        }
        json.append("]}");

        respond(ex, 200, "application/json", json.toString().getBytes());
    }

    // ── Helpers ────────────────────────────────────────────────────────────

    private void serveJson(HttpExchange ex, String path) throws IOException {
        File f = new File(path);
        if (!f.exists()) { respond(ex, 404, "application/json", "{\"error\":\"not found\"}".getBytes()); return; }
        respond(ex, 200, "application/json", Files.readAllBytes(f.toPath()));
    }

    private void respond(HttpExchange ex, int code, String mime, byte[] body) throws IOException {
        ex.getResponseHeaders().set("Content-Type", mime);
        ex.sendResponseHeaders(code, body.length);
        try (OutputStream os = ex.getResponseBody()) { os.write(body); }
    }

    private void cors(HttpExchange ex) {
        ex.getResponseHeaders().add("Access-Control-Allow-Origin",  "*");
        ex.getResponseHeaders().add("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
        ex.getResponseHeaders().add("Access-Control-Allow-Headers", "Content-Type");
    }

    private boolean preflight(HttpExchange ex) throws IOException {
        if ("OPTIONS".equalsIgnoreCase(ex.getRequestMethod())) {
            ex.sendResponseHeaders(204, -1); return true;
        }
        return false;
    }

    private static byte[] readAll(InputStream is) throws IOException {
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        byte[] tmp = new byte[8192]; int n;
        while ((n = is.read(tmp)) != -1) buf.write(tmp, 0, n);
        return buf.toByteArray();
    }

    private static String escape(String s) {
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }

    private static String extractStr(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return null;
        if (json.charAt(i) == '"') {
            StringBuilder sb = new StringBuilder(); i++;
            while (i < json.length()) {
                char c = json.charAt(i++);
                if (c == '"') break;
                if (c == '\\' && i < json.length()) { sb.append(json.charAt(i++)); continue; }
                sb.append(c);
            }
            return sb.toString();
        }
        if (json.startsWith("null", i)) return null;
        int end = i;
        while (end < json.length() && ",}]".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String guessMime(String path) {
        if (path.endsWith(".html")) return "text/html";
        if (path.endsWith(".js"))   return "application/javascript";
        if (path.endsWith(".json")) return "application/json";
        if (path.endsWith(".css"))  return "text/css";
        if (path.endsWith(".png"))  return "image/png";
        if (path.endsWith(".svg"))  return "image/svg+xml";
        return "application/octet-stream";
    }

    private static String folderPicker(String startDir) {
        final String[] result = {null};
        try {
            SwingUtilities.invokeAndWait(() -> {
                try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); } catch (Exception ignore) {}
                JFileChooser fc = new JFileChooser(startDir);
                fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                fc.setDialogTitle("Select output folder for PDF");
                if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
                    result[0] = fc.getSelectedFile().getAbsolutePath();
            });
        } catch (Exception e) {
            System.err.println("[LocalServer] Folder picker error: " + e.getMessage());
        }
        return result[0];
    }
}
