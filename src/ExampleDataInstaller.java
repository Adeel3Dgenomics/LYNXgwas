import java.io.*;
import java.net.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.util.*;
import java.util.zip.*;

/**
 * "Add example data": downloads the LYNXgwas article data (one zip, published on Zenodo), checks it, and
 * merges its projects and cross-dataset runs into the workspace (the current directory), with progress.
 *
 * Steps: find the latest version of the Zenodo record (or use a zip the user already has) -> check free
 * disk space -> download with resume -> verify size and checksum -> extract into example-data/staging,
 * refusing any entry that would land outside it -> fill in the workspace paths (pack.json lists the files
 * that hold placeholders) -> fetch the GRCh37 GENCODE annotation the projects use, if missing -> move the
 * projects and runs in. Anything already in the workspace with the same name is skipped, never replaced,
 * so the user's own projects are never touched.
 *
 * The pack is built by scripts/export_example_data.py. Its projects are marked results_only in their
 * project.json: everything can be browsed, but nothing that needs the original summary statistics can run.
 */
public class ExampleDataInstaller {

    /** Zenodo concept record of the example data; any version's record id also works. Overridable in CONFIG. */
    static final String DEFAULT_ZENODO_RECORD = "";
    static final String ZENODO_API = "https://zenodo.org/api/records/";
    static final File CONFIG = new File("config/example_data.json");
    static final File HOME = new File("example-data");
    static final File STATE = new File(HOME, "state.json");
    static final File STAGING = new File(HOME, "staging");
    static final String GFF3_NAME = "gencode.v37lift37.annotation.gff3";
    static final String GFF3_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/"
        + "GRCh37_mapping/gencode.v37lift37.annotation.gff3.gz";
    static final long GFF3_BYTES = 1_700_000_000L;
    /** The projects of the example-data pack (version 1), to tell whether a workspace already has them. */
    static final List<String> PACK_PROJECTS = Arrays.asList(
        "alz-demo", "alz-src-bellenguez2022", "alz-src-eadb2026", "alz-src-finngen-r12", "alz-src-nicolas2025",
        "cad-demo", "cad-src-aragam2022", "cad-src-hartiala2021", "cad-src-nam2026", "cad-src-nikpay2015",
        "ibd-demo", "ibd-src-delange2017uc", "ibd-src-finngen-r12", "ibd-src-liu2015", "ibd-src-zorina2023",
        "ra-demo", "ra-src-eyre2012", "ra-src-finngen-r12", "ra-src-glanville2021", "ra-src-stahl2010",
        "scz-demo", "scz-src-clozuk2018", "scz-src-pgc1-2011", "scz-src-pgc2-2014", "scz-src-sweden2013",
        "t2d-demo", "t2d-src-diamante2022", "t2d-src-finngen-r12", "t2d-src-mahajan2018", "t2d-src-scott2017");
    static final Set<String> ALLOWED_TOP =new HashSet<>(Arrays.asList(
        "projects", "output", "regulatory_data", "pack.json", "README.md", "LICENSE.txt"));

    /** Where the pack comes from: a Zenodo file or a configured URL. */
    public static class Source {
        String title = "LYNXgwas article data", version = "", url = "", fileName = "", md5 = "", sha256 = "",
               doi = "", page = "", error;
        long size;
    }

    public static class Progress {
        volatile String phase = "Starting", message = "", error;
        volatile long bytesDone, bytesTotal, startedAt = System.currentTimeMillis(), phaseStartedAt = startedAt;
        volatile int filesDone, filesTotal;
        volatile boolean done, cancelled;
        volatile String kind = "example";            // "example" (the article data) or "import" (any other pack)
        final List<String> installed = Collections.synchronizedList(new ArrayList<>());
        final List<String> skipped = Collections.synchronizedList(new ArrayList<>());

        void phase(String p, long total) {
            phase = p; bytesDone = 0; bytesTotal = total; filesDone = 0; filesTotal = 0;
            phaseStartedAt = System.currentTimeMillis();
        }

        String toJson() {
            long el = Math.max(1, System.currentTimeMillis() - phaseStartedAt);
            double speed = bytesDone * 1000.0 / el;
            StringBuilder j = new StringBuilder("{");
            j.append("\"phase\":\"").append(esc(phase)).append('"');
            j.append(",\"kind\":\"").append(esc(kind)).append('"');
            j.append(",\"message\":\"").append(esc(message)).append('"');
            j.append(",\"bytes_done\":").append(bytesDone).append(",\"bytes_total\":").append(bytesTotal);
            j.append(",\"files_done\":").append(filesDone).append(",\"files_total\":").append(filesTotal);
            j.append(",\"speed\":").append((long) speed);
            j.append(",\"done\":").append(done).append(",\"cancelled\":").append(cancelled);
            if (error != null) j.append(",\"error\":\"").append(esc(error)).append('"');
            j.append(",\"installed\":").append(list(installed)).append(",\"skipped\":").append(list(skipped));
            return j.append('}').toString();
        }
    }

    private volatile Progress job;
    private volatile Source cachedSource;
    private volatile long cachedAt;
    private final Runnable afterInstall;

    public ExampleDataInstaller(Runnable afterInstall) { this.afterInstall = afterInstall; }

    // ── Status, dismiss, cancel ─────────────────────────────────────────────

    /** {available, installed, dismissed, title, version, size, doi, job} for the home page. */
    public String statusJson() {
        Map<String, Object> st = readState();
        Source s = source(false);
        boolean running = job != null && !job.done;
        StringBuilder j = new StringBuilder("{");
        j.append("\"available\":").append(s.error == null);
        if (s.error != null) j.append(",\"source_error\":\"").append(esc(s.error)).append('"');
        j.append(",\"title\":\"").append(esc(s.title)).append('"');
        j.append(",\"version\":\"").append(esc(s.version)).append('"');
        j.append(",\"size\":").append(s.size);
        j.append(",\"doi\":\"").append(esc(s.doi)).append('"');
        j.append(",\"page\":\"").append(esc(s.page)).append('"');
        Object inst = st.get("installed_version");
        // The example projects may already be here without this installer (copied, or the workspace they came
        // from): when all of them are, there is nothing to offer
        int present = 0;
        for (String id : PACK_PROJECTS) if (new File("projects", id + "/config.properties").isFile()) present++;
        if (inst == null && present == PACK_PROJECTS.size()) inst = "present";
        j.append(",\"installed\":").append(inst == null ? "null" : "\"" + esc(String.valueOf(inst)) + "\"");
        j.append(",\"projects_present\":").append(present).append(",\"projects_total\":").append(PACK_PROJECTS.size());
        j.append(",\"dismissed\":").append(Boolean.TRUE.equals(st.get("dismissed")));
        j.append(",\"running\":").append(running);
        if (job != null) j.append(",\"job\":").append(job.toJson());
        return j.append('}').toString();
    }

    public String progressJson() { return job == null ? "{\"phase\":\"idle\",\"done\":true}" : job.toJson(); }

    public void cancel() { Progress p = job; if (p != null && !p.done) p.cancelled = true; }

    public synchronized void setDismissed(boolean dismissed) throws IOException {
        Map<String, Object> st = readState();
        st.put("dismissed", dismissed);
        writeState(st);
    }

    // ── Install ─────────────────────────────────────────────────────────────

    /** Starts the install in the background; localZip is a pack the user already has, or null to download. */
    public synchronized String start(String localZip) {
        if (job != null && !job.done) return "An install is already running.";
        if (localZip != null && !localZip.isEmpty() && !new File(localZip).isFile()) return "File not found: " + localZip;
        Progress p = new Progress();
        job = p;
        Thread t = new Thread(() -> run(p, localZip == null || localZip.isEmpty() ? null : new File(localZip)), "example-data");
        t.setDaemon(true);
        t.start();
        return null;
    }

    private void run(Progress p, File localZip) {
        File zip = null;
        boolean downloaded = false;
        try {
            HOME.mkdirs();
            Source s = null;
            if (localZip != null) {
                zip = localZip;
                p.kind = localZip.getName().startsWith("lynxgwas-example-data") ? "example" : "import";
            } else {
                p.phase("Finding the latest version", 0);
                s = source(true);
                if (s.error != null) throw new IOException(s.error);
                zip = new File(HOME, s.fileName);
                if (!(zip.isFile() && zip.length() == s.size)) {
                    long have = new File(HOME, s.fileName + ".part").length();
                    checkSpace(s.size - have + s.size * 3 + gffNeed(), "the download and its extraction");
                    download(p, s.url, new File(HOME, s.fileName + ".part"), s.size, "Downloading");
                    Files.move(new File(HOME, s.fileName + ".part").toPath(), zip.toPath(), StandardCopyOption.REPLACE_EXISTING);
                }
                downloaded = true;
                verify(p, zip, s);
            }

            Map<String, Object> pack;
            try (ZipFile zf = new ZipFile(zip)) {
                ZipEntry pe = zf.getEntry("pack.json");
                if (pe == null) throw new IOException("This file is not a LYNXgwas project pack (no pack.json).");
                try (InputStream in = zf.getInputStream(pe)) {
                    pack = MiniJson.asObject(MiniJson.parse(new String(readAll(in), StandardCharsets.UTF_8)));
                }
                if (!"1".equals(String.valueOf(num(pack.get("format")))))
                    throw new IOException("This pack needs a newer LYNXgwas (pack format " + pack.get("format") + ").");
                if ("projects".equals(MiniJson.getStr(pack, "kind", ""))) p.kind = "import";
                long unzipped = num(pack.get("uncompressed_bytes"));
                checkSpace(unzipped + gffNeed() + 500_000_000L, "extracting the pack");
                extract(p, zf, unzipped);
            }
            checkCancel(p);

            p.phase("Setting up the projects", 0);
            List<Object> rw = pack.get("rewrite_files") instanceof List ? MiniJson.asArray(pack.get("rewrite_files")) : new ArrayList<>();
            fillPaths(p, rw);
            checkCancel(p);

            ensureGff3(p, pack);
            checkCancel(p);

            p.phase("Adding the projects to LYNXgwas", 0);
            merge(p);
            restoreTimes(pack, p.installed);
            restamp(p.installed);

            if ("example".equals(p.kind)) {
                Map<String, Object> st = readState();
                st.put("installed_version", s != null ? s.version : String.valueOf(pack.getOrDefault("version", "")));
                st.put("installed_at", java.time.Instant.now().toString());
                if (s != null && !s.doi.isEmpty()) st.put("doi", s.doi);
                st.put("projects", new ArrayList<>(p.installed));
                st.put("dismissed", true);
                writeState(st);
            }

            deleteUnder(STAGING);
            if (downloaded && zip != null) zip.delete();
            if (afterInstall != null) afterInstall.run();
            p.phase("Done", 0);
            p.message = p.installed.size() + " added" + (p.skipped.isEmpty() ? "" : ", " + p.skipped.size() + " already present (kept as they were)");
        } catch (CancelledException e) {
            p.phase("Cancelled", 0);
            p.cancelled = true;
            p.message = "Cancelled. A partial download is kept, so the next try continues where this one stopped.";
            try { deleteUnder(STAGING); } catch (IOException ignored) {}
        } catch (Exception e) {
            p.error = e.getMessage() == null ? e.toString() : e.getMessage();
            p.phase("Failed", 0);
            try { deleteUnder(STAGING); } catch (IOException ignored) {}
            System.err.println("[Import] failed: " + p.error);
        } finally {
            p.done = true;
        }
    }

    // ── Source (Zenodo) ─────────────────────────────────────────────────────

    private Source source(boolean fresh) {
        Source c = cachedSource;
        if (!fresh && c != null && System.currentTimeMillis() - cachedAt < 10 * 60_000L) return c;
        Source s = new Source();
        try {
            Map<String, Object> cfg = CONFIG.isFile()
                ? MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(CONFIG.toPath()), StandardCharsets.UTF_8)))
                : new HashMap<>();
            String record = MiniJson.getStr(cfg, "zenodo_record", DEFAULT_ZENODO_RECORD).trim();
            if (!record.isEmpty()) {
                fromZenodo(s, MiniJson.getStr(cfg, "zenodo_api", ZENODO_API) + record + "/versions/latest");
            } else if (!MiniJson.getStr(cfg, "url", "").isEmpty()) {
                s.url = MiniJson.getStr(cfg, "url", "");
                s.fileName = MiniJson.getStr(cfg, "file", s.url.replaceAll(".*/", ""));
                s.size = num(cfg.get("size"));
                s.md5 = MiniJson.getStr(cfg, "md5", "");
                s.sha256 = MiniJson.getStr(cfg, "sha256", "");
                s.version = MiniJson.getStr(cfg, "version", "");
                s.title = MiniJson.getStr(cfg, "title", s.title);
                s.doi = MiniJson.getStr(cfg, "doi", "");
            } else {
                s.error = "The example data is not published yet. Use \"Install from a file\" with a pack you already have.";
            }
            if (s.error == null && (s.url.isEmpty() || s.size <= 0))
                s.error = "The example-data record does not list a pack file.";
            if (s.error == null && !s.fileName.matches("[A-Za-z0-9._-]+\\.zip"))
                s.error = "Unexpected pack file name: " + s.fileName;
        } catch (UnknownHostException | ConnectException | SocketTimeoutException e) {
            s.error = "Could not reach Zenodo to find the example data (" + e.getMessage() + "). Check the internet connection.";
        } catch (Exception e) {
            s.error = "Could not read the example-data record on Zenodo: " + e.getMessage();
        }
        cachedSource = s;
        cachedAt = System.currentTimeMillis();
        return s;
    }

    /** Reads a Zenodo record (latest version) and picks its lynxgwas-example-data-*.zip file. */
    private static void fromZenodo(Source s, String api) throws IOException {
        HttpURLConnection c = open(api, 0, "application/json");
        if (c.getResponseCode() != 200) throw new IOException("Zenodo answered " + c.getResponseCode());
        Map<String, Object> rec;
        try (InputStream in = c.getInputStream()) { rec = MiniJson.asObject(MiniJson.parse(new String(readAll(in), StandardCharsets.UTF_8))); }
        Map<String, Object> meta = rec.get("metadata") instanceof Map ? MiniJson.asObject(rec.get("metadata")) : new HashMap<>();
        s.title = MiniJson.getStr(meta, "title", s.title);
        s.version = MiniJson.getStr(meta, "version", "");
        s.doi = MiniJson.getStr(rec, "doi", "");
        Map<String, Object> links = rec.get("links") instanceof Map ? MiniJson.asObject(rec.get("links")) : new HashMap<>();
        s.page = MiniJson.getStr(links, "self_html", MiniJson.getStr(links, "html", ""));
        List<Object> files = new ArrayList<>();
        Object f = rec.get("files");
        if (f instanceof List) files = MiniJson.asArray(f);
        else if (f instanceof Map && MiniJson.asObject(f).get("entries") instanceof Map)
            files = new ArrayList<>(MiniJson.asObject(MiniJson.asObject(f).get("entries")).values());
        for (Object o : files) {
            Map<String, Object> e = MiniJson.asObject(o);
            String key = MiniJson.getStr(e, "key", MiniJson.getStr(e, "filename", ""));
            if (!key.matches("lynxgwas-example-data.*\\.zip")) continue;
            s.fileName = key;
            s.size = num(e.containsKey("size") ? e.get("size") : e.get("filesize"));
            String sum = MiniJson.getStr(e, "checksum", "");
            if (sum.startsWith("md5:")) s.md5 = sum.substring(4);
            Map<String, Object> l = e.get("links") instanceof Map ? MiniJson.asObject(e.get("links")) : new HashMap<>();
            s.url = MiniJson.getStr(l, "content", MiniJson.getStr(l, "self", MiniJson.getStr(l, "download", "")));
            if (!s.url.isEmpty() && !s.url.endsWith("/content") && s.url.contains("/api/records/")) s.url += "/content";
            return;
        }
        throw new IOException("the Zenodo record has no lynxgwas-example-data zip");
    }

    // ── Download, verify, extract ───────────────────────────────────────────

    /** Downloads url into part, continuing an earlier partial download when the server allows it. */
    private static void download(Progress p, String url, File part, long expected, String label) throws IOException {
        long have = part.isFile() ? part.length() : 0;
        if (expected > 0 && have > expected) { part.delete(); have = 0; }
        p.phase(label, expected);
        p.bytesDone = have;
        for (int attempt = 1; ; attempt++) {
            try {
                HttpURLConnection c = open(url, have);
                int code = c.getResponseCode();
                if (code == 416 && expected > 0 && have == expected) return;     // already complete
                boolean append = code == 206;
                if (code != 200 && code != 206) throw new IOException("the server answered " + code + " for " + url);
                if (!append) { have = 0; p.bytesDone = 0; }
                long len = c.getContentLengthLong();
                if (expected <= 0) { p.bytesTotal = len + have; }
                else if (len > 0 && have + len != expected) {
                    // retrying cannot help: the server holds a different file than the record describes
                    part.delete();
                    throw new MismatchException("The file on the server (" + (have + len) + " bytes) is not the one the "
                        + "example-data record describes (" + expected + " bytes). Please try again later.");
                }
                try (InputStream in = new BufferedInputStream(c.getInputStream(), 1 << 16);
                     OutputStream out = new BufferedOutputStream(new FileOutputStream(part, append), 1 << 16)) {
                    byte[] buf = new byte[1 << 16];
                    int n;
                    while ((n = in.read(buf)) > 0) {
                        out.write(buf, 0, n);
                        have += n;
                        p.bytesDone = have;
                        if (p.cancelled) throw new CancelledException();
                    }
                }
                if (expected > 0 && have != expected)
                    throw new IOException("the download stopped at " + have + " of " + expected + " bytes");
                return;
            } catch (CancelledException | MismatchException e) {
                throw e;
            } catch (IOException e) {
                if (attempt >= 5) throw new IOException("Download failed after 5 tries: " + e.getMessage());
                p.message = "Connection problem (" + e.getMessage() + "); retrying (" + attempt + "/5)";
                try { Thread.sleep(3000L * attempt); } catch (InterruptedException ie) { throw new CancelledException(); }
                have = part.isFile() ? part.length() : 0;
                p.message = "";
            }
        }
    }

    private static HttpURLConnection open(String url, long from) throws IOException { return open(url, from, null); }

    /** Opens url (following up to 6 HTTPS redirects), asking for bytes from `from` on when it is above 0. */
    private static HttpURLConnection open(String url, long from, String accept) throws IOException {
        String u = url;
        for (int hop = 0; hop < 6; hop++) {
            HttpURLConnection c = (HttpURLConnection) new URL(u).openConnection();
            c.setInstanceFollowRedirects(false);
            c.setConnectTimeout(20_000);
            c.setReadTimeout(120_000);
            c.setRequestProperty("User-Agent", "LYNXgwas/" + Config.PIPELINE_VERSION + " (example-data installer)");
            if (from > 0) c.setRequestProperty("Range", "bytes=" + from + "-");
            if (accept != null) c.setRequestProperty("Accept", accept);
            int code = c.getResponseCode();
            if (code >= 300 && code < 400 && c.getHeaderField("Location") != null) {
                u = new URL(new URL(u), c.getHeaderField("Location")).toString();
                if (!u.startsWith("https://")) throw new IOException("refusing a non-HTTPS redirect: " + u);
                c.disconnect();
                continue;
            }
            return c;
        }
        throw new IOException("too many redirects for " + url);
    }

    private static void verify(Progress p, File zip, Source s) throws IOException {
        if (s.size > 0 && zip.length() != s.size) {
            zip.delete();
            throw new IOException("The downloaded file has the wrong size; it was deleted. Please try again.");
        }
        if (s.md5.isEmpty() && s.sha256.isEmpty()) return;
        p.phase("Checking the download", zip.length());
        try {
            MessageDigest md5 = MessageDigest.getInstance("MD5"), sha = MessageDigest.getInstance("SHA-256");
            try (InputStream in = new BufferedInputStream(new FileInputStream(zip), 1 << 16)) {
                byte[] buf = new byte[1 << 16];
                int n;
                while ((n = in.read(buf)) > 0) {
                    md5.update(buf, 0, n); sha.update(buf, 0, n);
                    p.bytesDone += n;
                    if (p.cancelled) throw new CancelledException();
                }
            }
            if ((!s.md5.isEmpty() && !s.md5.equalsIgnoreCase(hex(md5.digest())))
                    || (!s.sha256.isEmpty() && !s.sha256.equalsIgnoreCase(hex(sha.digest())))) {
                zip.delete();
                throw new IOException("The downloaded file is damaged (checksum mismatch); it was deleted. Please try again.");
            }
        } catch (java.security.NoSuchAlgorithmException e) {
            throw new IOException(e);
        }
    }

    /** Extracts into STAGING, refusing entries outside it or outside the pack's known folders. */
    private static void extract(Progress p, ZipFile zf, long declared) throws IOException {
        deleteUnder(STAGING);
        Path root = STAGING.toPath().toAbsolutePath().normalize();
        Files.createDirectories(root);
        long limit = declared > 0 ? declared + declared / 10 + 50_000_000L : Long.MAX_VALUE, written = 0;
        p.phase("Unpacking", declared);
        p.filesTotal = zf.size();
        byte[] buf = new byte[1 << 16];
        for (Enumeration<? extends ZipEntry> en = zf.entries(); en.hasMoreElements(); ) {
            ZipEntry e = en.nextElement();
            String name = e.getName().replace('\\', '/');
            String top = name.contains("/") ? name.substring(0, name.indexOf('/')) : name;
            if (name.startsWith("/") || name.contains(":") || !ALLOWED_TOP.contains(top))
                throw new IOException("Refusing unexpected entry in the pack: " + name);
            Path out = root.resolve(name).normalize();
            if (!out.startsWith(root)) throw new IOException("Refusing entry outside the pack folder: " + name);
            if (e.isDirectory()) { Files.createDirectories(out); continue; }
            Files.createDirectories(out.getParent());
            try (InputStream in = zf.getInputStream(e); OutputStream os = new BufferedOutputStream(Files.newOutputStream(out), 1 << 16)) {
                int n;
                while ((n = in.read(buf)) > 0) {
                    os.write(buf, 0, n);
                    written += n;
                    p.bytesDone = written;
                    if (written > limit) throw new IOException("The pack unpacks to more than it declares; stopped.");
                }
            }
            if (e.getTime() > 0) out.toFile().setLastModified(e.getTime());
            p.filesDone++;
            if (p.cancelled) throw new CancelledException();
        }
    }

    /** Replaces the pack's path placeholders with this workspace's folders. */
    private static void fillPaths(Progress p, List<Object> files) throws IOException {
        String ws = new File(".").getCanonicalPath().replace('\\', '/');
        String[][] map = {
            {"@@LYNX_WORKSPACE@@", ws}, {"@@LYNX_APP@@", ws}, {"@@LYNX_RESOURCES@@", ws + "/resources"},
            {"@@LYNX_REF_PANEL@@", ws + "/refpanels"}, {"@@LYNX_GWAS@@", ws + "/example-data/gwas"}};
        Path root = STAGING.toPath().toAbsolutePath().normalize();
        p.filesTotal = files.size();
        for (Object o : files) {
            Path f = root.resolve(String.valueOf(o)).normalize();
            if (!f.startsWith(root) || !Files.isRegularFile(f)) { p.filesDone++; continue; }
            String t = new String(Files.readAllBytes(f), StandardCharsets.UTF_8), n = t;
            for (String[] m : map) n = n.replace(m[0], m[1]);
            if (!n.equals(t)) Files.write(f, n.getBytes(StandardCharsets.UTF_8));
            p.filesDone++;
            if (p.cancelled) throw new CancelledException();
        }
    }

    private static long gffNeed() { return new File("resources", GFF3_NAME).isFile() ? 0 : GFF3_BYTES + 60_000_000L; }

    /**
     * Fetches from GENCODE each gene annotation the pack's projects use that this workspace lacks
     * (pack.json "gff3": [{file, url}]; packs without it are the example data, which uses v37lift37).
     */
    private static void ensureGff3(Progress p, Map<String, Object> pack) throws IOException {
        List<String[]> wanted = new ArrayList<>();
        Object g = pack.get("gff3");
        if (g instanceof List) {
            for (Object o : MiniJson.asArray(g)) {
                Map<String, Object> m = MiniJson.asObject(o);
                wanted.add(new String[]{MiniJson.getStr(m, "file", ""), MiniJson.getStr(m, "url", "")});
            }
        } else if (g instanceof Map) {
            wanted.add(new String[]{MiniJson.getStr(MiniJson.asObject(g), "file", GFF3_NAME), MiniJson.getStr(MiniJson.asObject(g), "url", GFF3_URL)});
        } else {
            wanted.add(new String[]{GFF3_NAME, GFF3_URL});
        }
        for (String[] w : wanted) {
            if (!w[0].matches("[A-Za-z0-9._-]+\\.gff3")) continue;
            if (new File("resources", w[0]).isFile()) continue;
            if (w[1].isEmpty() || !w[1].startsWith("https://")) {
                p.message = "The gene annotation " + w[0] + " is not available to download; put it in resources/ (Edit to point a project elsewhere).";
                continue;
            }
            fetchGff3(p, w[0], w[1]);
        }
    }

    private static void fetchGff3(Progress p, String name, String url) throws IOException {
        File gff = new File("resources", name);
        new File("resources").mkdirs();
        File gz = new File(HOME, name + ".gz.part");
        download(p, url, gz, -1, "Downloading the gene annotation (GENCODE)");
        p.phase("Unpacking the gene annotation", gz.length());
        File tmp = new File("resources", name + ".part");
        try (CountingInputStream cin = new CountingInputStream(new FileInputStream(gz), p);
             InputStream in = new GZIPInputStream(new BufferedInputStream(cin, 1 << 16), 1 << 16);
             OutputStream out = new BufferedOutputStream(new FileOutputStream(tmp), 1 << 16)) {
            byte[] buf = new byte[1 << 16];
            int n;
            while ((n = in.read(buf)) > 0) {
                out.write(buf, 0, n);
                if (p.cancelled) throw new CancelledException();
            }
        }
        Files.move(tmp.toPath(), gff.toPath(), StandardCopyOption.REPLACE_EXISTING);
        gz.delete();
    }

    /** Moves the staged projects and runs into the workspace; anything with the same name is left alone. */
    private static void merge(Progress p) throws IOException {
        Path st = STAGING.toPath();
        String[][] groups = {{"projects", "projects"}, {"output/multi_locus", "output/multi_locus"}, {"output/mvsusie", "output/mvsusie"}};
        for (String[] g : groups) {
            File src = st.resolve(g[0]).toFile();
            File[] items = src.listFiles(File::isDirectory);
            if (items == null) continue;
            new File(g[1]).mkdirs();
            for (File it : items) {
                File dest = new File(g[1], it.getName());
                boolean isProject = g[0].equals("projects");
                if (dest.exists()) {
                    if (isProject) p.skipped.add(it.getName());
                    continue;
                }
                move(it.toPath(), dest.toPath());
                if (isProject) p.installed.add(it.getName());
            }
        }
        if (!p.installed.isEmpty()) removePipPlaceholder();
        File reg = st.resolve("regulatory_data").toFile();
        if (reg.isDirectory()) {
            Path from = reg.toPath();
            try (java.util.stream.Stream<Path> w = Files.walk(from)) {
                for (Path f : (Iterable<Path>) w.filter(Files::isRegularFile)::iterator) {
                    Path dest = Paths.get("regulatory_data").resolve(from.relativize(f).toString());
                    if (Files.exists(dest)) continue;
                    Files.createDirectories(dest.getParent());
                    move(f, dest);
                }
            }
        }
    }

    /**
     * The pip launcher seeds an empty projects/getting-started (the app will not start with no project at
     * all). Once real projects exist it only shows as a broken row, so it is removed, but only while it is
     * still exactly what the launcher wrote: a config and a header-only placeholder file, nothing else.
     */
    private static void removePipPlaceholder() {
        File dir = new File("projects", "getting-started");
        File[] files = dir.listFiles();
        if (files == null || files.length != 2) return;
        File cfg = new File(dir, "config.properties"), gwas = new File(dir, "placeholder_gwas.tsv");
        if (!cfg.isFile() || !gwas.isFile() || gwas.length() > 100) return;
        try {
            if (!new String(Files.readAllBytes(cfg.toPath()), StandardCharsets.UTF_8).contains("placeholder_gwas.tsv")) return;
            Files.delete(gwas.toPath());
            Files.delete(cfg.toPath());
            Files.delete(dir.toPath());
        } catch (IOException ignored) {}
    }

    /**
     * A project that arrived with its summary statistics was up to date where it was exported, but its
     * fingerprint names that machine's paths: stamp it again here so it is not shown as needing reprocessing.
     * Results-only projects are never stale and need nothing.
     */
    private static void restamp(List<String> ids) {
        synchronized (ids) {
            for (String id : ids) {
                String dir = new File("projects", id).getAbsolutePath();
                try {
                    ProjectMetadata pm = ProjectMetadata.load(dir);
                    if (pm == null || pm.resultsOnly) continue;
                    Config cfg = Config.loadFromProject(dir);
                    pm.coreInputFingerprint = ProjectMetadata.computeCoreInputFingerprint(cfg, dir);
                    pm.annotationFingerprint = ProjectMetadata.computeAnnotationFingerprint(new File(dir, "annotations.yaml").getAbsolutePath());
                    String last = pm.lastProcessed;
                    pm.save(dir);
                    if (last != null && !last.isEmpty()) {          // save() stamps "now"; keep when it was processed
                        File f = new File(dir, "project.json");
                        String t = new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
                        Files.write(f.toPath(), t.replaceFirst("\"last_processed\": \"[^\"]*\"",
                            java.util.regex.Matcher.quoteReplacement("\"last_processed\": \"" + last + "\"")).getBytes(StandardCharsets.UTF_8));
                    }
                } catch (Exception e) {
                    System.err.println("[Import] could not re-stamp " + id + ": " + e.getMessage());
                }
            }
        }
    }

    /** Gives files their exact modification time from pack.json (a zip keeps it only to 2 s). */
    private static void restoreTimes(Map<String, Object> pack, List<String> installed) {
        if (!(pack.get("mtimes") instanceof Map)) return;
        for (Map.Entry<String, Object> e : MiniJson.asObject(pack.get("mtimes")).entrySet()) {
            String arc = e.getKey();
            String[] parts = arc.split("/");
            if (parts.length < 3 || !parts[0].equals("projects") || !installed.contains(parts[1]) || arc.contains("..")) continue;
            File f = new File(arc);
            long t = num(e.getValue());
            if (f.isFile() && t > 0) f.setLastModified(t);
        }
    }

    private static void move(Path from, Path to) throws IOException {
        try {
            Files.move(from, to, StandardCopyOption.ATOMIC_MOVE);
        } catch (IOException e) {
            Files.move(from, to);
        }
    }

    // ── Helpers ─────────────────────────────────────────────────────────────

    private static void checkSpace(long need, String what) throws IOException {
        long free = new File(".").getAbsoluteFile().getUsableSpace();
        if (need > 0 && free < need)
            throw new IOException(String.format(Locale.ROOT, "Not enough disk space for %s: %.1f GB needed, %.1f GB free.",
                what, need / 1e9, free / 1e9));
    }

    private static void checkCancel(Progress p) { if (p.cancelled) throw new CancelledException(); }

    /** Deletes a folder, but only inside example-data/. */
    private static void deleteUnder(File dir) throws IOException {
        if (!dir.exists()) return;
        Path home = HOME.toPath().toAbsolutePath().normalize(), d = dir.toPath().toAbsolutePath().normalize();
        if (!d.startsWith(home) || d.equals(home)) throw new IOException("refusing to delete " + d);
        try (java.util.stream.Stream<Path> w = Files.walk(d)) {
            List<Path> all = new ArrayList<>();
            w.forEach(all::add);
            Collections.reverse(all);
            for (Path x : all) Files.deleteIfExists(x);
        }
    }

    private static Map<String, Object> readState() {
        try {
            if (STATE.isFile())
                return MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(STATE.toPath()), StandardCharsets.UTF_8)));
        } catch (Exception ignored) {}
        return new LinkedHashMap<>();
    }

    private static void writeState(Map<String, Object> st) throws IOException {
        HOME.mkdirs();
        Files.write(STATE.toPath(), MiniJson.encode(st).getBytes(StandardCharsets.UTF_8));
    }

    static long num(Object o) {
        if (o instanceof Number) return ((Number) o).longValue();
        try { return o == null ? 0 : Long.parseLong(String.valueOf(o).trim()); } catch (NumberFormatException e) { return 0; }
    }

    private static byte[] readAll(InputStream in) throws IOException {
        ByteArrayOutputStream b = new ByteArrayOutputStream();
        byte[] buf = new byte[8192];
        int n;
        while ((n = in.read(buf)) > 0) b.write(buf, 0, n);
        return b.toByteArray();
    }

    private static String hex(byte[] d) {
        StringBuilder sb = new StringBuilder();
        for (byte x : d) sb.append(String.format("%02x", x & 0xff));
        return sb.toString();
    }

    private static String list(List<String> l) {
        StringBuilder j = new StringBuilder("[");
        synchronized (l) {
            for (int i = 0; i < l.size(); i++) j.append(i > 0 ? "," : "").append('"').append(esc(l.get(i))).append('"');
        }
        return j.append(']').toString();
    }

    static String esc(String s) {
        if (s == null) return "";
        StringBuilder b = new StringBuilder();
        for (char c : s.toCharArray()) {
            if (c == '"' || c == '\\') b.append('\\').append(c);
            else if (c < 0x20) b.append(String.format("\\u%04x", (int) c));
            else b.append(c);
        }
        return b.toString();
    }

    static class CancelledException extends RuntimeException {}

    /** A download that cannot succeed by retrying. */
    static class MismatchException extends IOException { MismatchException(String m) { super(m); } }

    /** Counts compressed bytes read, for the gunzip progress bar. */
    static class CountingInputStream extends FilterInputStream {
        final Progress p;
        CountingInputStream(InputStream in, Progress p) { super(in); this.p = p; }
        @Override public int read() throws IOException { int r = super.read(); if (r >= 0) p.bytesDone++; return r; }
        @Override public int read(byte[] b, int off, int len) throws IOException {
            int n = super.read(b, off, len);
            if (n > 0) p.bytesDone += n;
            return n;
        }
    }
}
