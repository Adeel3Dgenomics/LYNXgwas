import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.zip.*;

/**
 * Export: packs one or more projects into a zip that "Import" (ExampleDataInstaller) reads back, on this
 * machine or another. Same format as the example data (scripts/export_example_data.py): the latest run
 * of each tool at each locus without its bulky inputs, the viewer data, heritability and the other saved
 * results; optionally the cross-dataset runs whose datasets are all in the export (with their joint
 * fine-mapping runs) and the projects' GWAS summary statistics. Without the summary statistics the
 * projects arrive as results only. Machine paths become placeholders the importer fills in.
 *
 * The zip is written to exports/ in the workspace and offered for download.
 */
public class ProjectPackExporter {

    static final File OUT_DIR = new File("exports");
    static final Set<String> TEXT_EXT = new HashSet<>(Arrays.asList(
        "json", "properties", "tsv", "txt", "log", "r", "yaml", "yml", "js", "csv", "md", "key", "snps", "cojo", "master", "html"));
    static final java.util.regex.Pattern RUN_EXCLUDE = java.util.regex.Pattern.compile(
        "^(finemap\\.ld|finemap\\.z|ref_bigsnp\\..*|susie_input\\.tsv|cojo_input\\.ma|ma_audit\\.tsv|magma_snploc\\.txt|"
        + "magma_input\\.pval|magma_step\\d+\\.genes\\.annot\\.txt|.*\\.cma\\.cojo|.*\\.(bed|bim|fam|bk|rds|ld))$",
        java.util.regex.Pattern.CASE_INSENSITIVE);
    static final Set<String> MVSUSIE_KEEP = new HashSet<>(Arrays.asList(
        "meta.json", "result.json", "manifest.json", "run.log", "extract.txt", "ukb_snps.tsv"));
    static final Set<String> PROJECT_TOP_SKIP = new HashSet<>(Arrays.asList(
        ".fingerprint_cache", "plink_subsets", "ld_results", "tmp", "plots", "loci_analysis", "input"));
    /** Gene annotations the importer can fetch from GENCODE when the importing workspace lacks them. */
    static final Map<String, String> KNOWN_GFF3 = new LinkedHashMap<>();
    static {
        KNOWN_GFF3.put("gencode.v37lift37.annotation.gff3", ExampleDataInstaller.GFF3_URL);
        KNOWN_GFF3.put("gencode.v44.basic.annotation.gff3",
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gff3.gz");
    }

    private volatile ExampleDataInstaller.Progress job;
    private volatile String lastFile;

    public String progressJson() {
        ExampleDataInstaller.Progress p = job;
        if (p == null) return "{\"phase\":\"idle\",\"done\":true}";
        String j = p.toJson();
        return lastFile != null && p.done && p.error == null
            ? j.substring(0, j.length() - 1) + ",\"file\":\"" + ExampleDataInstaller.esc(lastFile) + "\",\"folder\":\""
                + ExampleDataInstaller.esc(OUT_DIR.getAbsolutePath().replace('\\', '/')) + "\"}"
            : j;
    }

    public void cancel() { ExampleDataInstaller.Progress p = job; if (p != null && !p.done) p.cancelled = true; }

    /** The finished export with this name, or null (names are checked, so nothing outside exports/ is served). */
    public static File exported(String name) {
        if (name == null || !name.matches("lynxgwas-[A-Za-z0-9._-]+\\.zip")) return null;
        File f = new File(OUT_DIR, name);
        return f.isFile() ? f : null;
    }

    public synchronized String start(List<String> ids, boolean includeGwas, boolean includeRuns) {
        if (job != null && !job.done) return "An export is already running.";
        if (ids == null || ids.isEmpty()) return "Choose at least one project.";
        for (String id : ids)
            if (!id.matches("[A-Za-z0-9._-]+") || !new File("projects", id + "/config.properties").isFile())
                return "Project not found: " + id;
        ExampleDataInstaller.Progress p = new ExampleDataInstaller.Progress();
        job = p;
        lastFile = null;
        Thread t = new Thread(() -> run(p, new ArrayList<>(ids), includeGwas, includeRuns), "project-export");
        t.setDaemon(true);
        t.start();
        return null;
    }

    private void run(ExampleDataInstaller.Progress p, List<String> ids, boolean includeGwas, boolean includeRuns) {
        File part = null;
        try {
            p.phase("Collecting the files", 0);
            String ws = new File(".").getCanonicalPath().replace('\\', '/');
            List<String[]> rules = new ArrayList<>();                    // {absolute path spelling, placeholder}
            addRule(rules, ws, "@@LYNX_WORKSPACE@@");
            Map<File, String> files = new LinkedHashMap<>();              // source -> path in the pack
            Set<String> gffNames = new LinkedHashSet<>();
            Map<String, Object> mtimes = new LinkedHashMap<>();           // exact times a zip cannot hold (2 s)
            Set<String> eids = new LinkedHashSet<>();
            boolean allGwas = true;

            for (String id : ids) {
                Config cfg = Config.loadFromProject(new File("projects", id).getAbsolutePath());
                File gwas = new File(cfg.gwasFile);
                boolean withGwas = includeGwas && gwas.isFile();
                allGwas &= withGwas;
                if (withGwas) {
                    files.put(gwas, "projects/" + id + "/input/" + gwas.getName());
                    // saved checks (sample size, enrichment, QC) are keyed on the file's size and time
                    mtimes.put("projects/" + id + "/input/" + gwas.getName(), gwas.lastModified());
                    addRule(rules, gwas.getAbsolutePath(), "@@LYNX_WORKSPACE@@/projects/" + id + "/input/" + gwas.getName());
                }
                if (gwas.getParentFile() != null) addRule(rules, gwas.getAbsoluteFile().getParent(), "@@LYNX_GWAS@@");
                if (!cfg.refPanelPath.isEmpty() && new File(cfg.refPanelPath).getAbsoluteFile().getParent() != null)
                    addRule(rules, new File(cfg.refPanelPath).getAbsoluteFile().getParent(), "@@LYNX_REF_PANEL@@");
                if (!cfg.gff3File.isEmpty()) {
                    File g = new File(cfg.gff3File).getAbsoluteFile();
                    gffNames.add(g.getName());
                    if (g.getParent() != null) addRule(rules, g.getParent(), "@@LYNX_RESOURCES@@");
                }
                collectProject(new File("projects", id), "projects/" + id, files);
                String eid = RegulatoryPeakIndex.resolveEid(id);
                if (eid != null) eids.add(eid);
                checkCancel(p);
            }
            rules.sort((a, b) -> b[0].length() - a[0].length());

            List<String> jobs = new ArrayList<>();
            if (includeRuns) collectRuns(new HashSet<>(ids), files, jobs);
            // the Roadmap histone peaks of the projects' tissues: regulatory tracks and enrichment need them
            for (String eid : eids)
                if (eid.matches("E\\d+")) addTree(new File("regulatory_data", eid), "regulatory_data/" + eid, files, f -> true);

            long total = 0;
            for (File f : files.keySet()) total += f.length();
            OUT_DIR.mkdirs();
            String stamp = new SimpleDateFormat("yyyyMMdd-HHmm").format(new Date());
            String name = ids.size() == 1 ? "lynxgwas-project-" + ids.get(0) + "-" + stamp + ".zip"
                                          : "lynxgwas-projects-" + ids.size() + "-" + stamp + ".zip";
            long free = OUT_DIR.getAbsoluteFile().getUsableSpace();
            if (free < total / 2 + 200_000_000L)
                throw new IOException(String.format(Locale.ROOT, "Not enough disk space for the export: about %.1f GB needed, %.1f GB free.",
                    (total / 2 + 200_000_000L) / 1e9, free / 1e9));

            part = new File(OUT_DIR, name + ".part");
            List<String> rewritten = new ArrayList<>();
            p.phase("Writing " + name, total);
            p.filesTotal = files.size();
            try (ZipOutputStream z = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(part), 1 << 16))) {
                z.setLevel(6);
                for (Map.Entry<File, String> e : files.entrySet()) {
                    File src = e.getKey();
                    String arc = e.getValue();
                    ZipEntry ze = new ZipEntry(arc);
                    ze.setTime(src.lastModified());
                    z.putNextEntry(ze);
                    if (isText(arc)) {
                        String t = new String(Files.readAllBytes(src.toPath()), StandardCharsets.UTF_8);
                        for (String[] r : rules) if (t.contains(r[0])) t = t.replace(r[0], r[1]);
                        t = projectOverride(arc, t, ids, includeGwas);
                        if (t.contains("@@LYNX_")) rewritten.add(arc);
                        z.write(t.getBytes(StandardCharsets.UTF_8));
                        p.bytesDone += src.length();
                    } else {
                        try (InputStream in = new BufferedInputStream(new FileInputStream(src), 1 << 16)) {
                            byte[] buf = new byte[1 << 16];
                            int n;
                            while ((n = in.read(buf)) > 0) { z.write(buf, 0, n); p.bytesDone += n; checkCancel(p); }
                        }
                    }
                    z.closeEntry();
                    p.filesDone++;
                    checkCancel(p);
                }
                Map<String, Object> pack = new LinkedHashMap<>();
                pack.put("format", 1);
                pack.put("name", name.replaceAll("\\.zip$", ""));
                pack.put("kind", "projects");
                pack.put("created", java.time.Instant.now().toString());
                pack.put("projects", new ArrayList<Object>(ids));
                pack.put("cross_dataset_runs", new ArrayList<Object>(jobs));
                pack.put("includes_gwas", includeGwas && allGwas);
                pack.put("files", files.size());
                pack.put("uncompressed_bytes", total);
                pack.put("rewrite_files", new ArrayList<Object>(rewritten));
                pack.put("mtimes", mtimes);
                List<Object> gff = new ArrayList<>();
                for (String g : gffNames) {
                    Map<String, Object> m = new LinkedHashMap<>();
                    m.put("file", g);
                    if (KNOWN_GFF3.containsKey(g)) m.put("url", KNOWN_GFF3.get(g));
                    gff.add(m);
                }
                pack.put("gff3", gff);
                z.putNextEntry(new ZipEntry("pack.json"));
                z.write(MiniJson.encode(pack).getBytes(StandardCharsets.UTF_8));
                z.closeEntry();
            }
            File out = new File(OUT_DIR, name);
            Files.move(part.toPath(), out.toPath(), StandardCopyOption.REPLACE_EXISTING);
            lastFile = name;
            p.phase("Done", 0);
            p.message = String.format(Locale.ROOT, "%d project%s%s, %.0f MB", ids.size(), ids.size() == 1 ? "" : "s",
                jobs.isEmpty() ? "" : " and " + jobs.size() + " cross-dataset run" + (jobs.size() == 1 ? "" : "s"), out.length() / 1e6);
        } catch (ExampleDataInstaller.CancelledException e) {
            p.phase("Cancelled", 0);
            p.cancelled = true;
            p.message = "Export cancelled.";
            if (part != null) part.delete();
        } catch (Exception e) {
            p.error = e.getMessage() == null ? e.toString() : e.getMessage();
            p.phase("Failed", 0);
            if (part != null) part.delete();
        } finally {
            p.done = true;
        }
    }

    /** What a project contributes: everything the app reads, with the latest run of each tool at each locus. */
    private static void collectProject(File dir, String arc, Map<File, String> out) throws IOException {
        File[] top = dir.listFiles();
        if (top == null) return;
        for (File e : top) {
            String n = e.getName();
            if (PROJECT_TOP_SKIP.contains(n) || n.endsWith(".bak") || n.endsWith(".hg18") || n.endsWith("_done")) continue;
            if (e.isFile()) out.put(e, arc + "/" + n);
            else if (e.isDirectory()) addTree(e, arc + "/" + n, out, n.equals("data") ? f -> !f.matches("locus_\\d+\\.js") : f -> true);
        }
        File[] loci = new File(dir, "loci_analysis").listFiles(File::isDirectory);
        if (loci == null) return;
        for (File L : loci) {
            String la = arc + "/loci_analysis/" + L.getName();
            if (new File(L, "build.log").isFile()) out.put(new File(L, "build.log"), la + "/build.log");
            for (String n : new String[]{"consistency_summary.json", "consistency_report.tsv"})
                if (new File(L, "ld/" + n).isFile()) out.put(new File(L, "ld/" + n), la + "/ld/" + n);
            if (new File(L, "failures").isDirectory()) addTree(new File(L, "failures"), la + "/failures", out, f -> true);
            for (File run : latestRuns(new File(L, "runs")))
                addTree(run, la + "/runs/" + run.getName(), out, f -> !RUN_EXCLUDE.matcher(f).matches());
        }
    }

    /** Cross-dataset runs whose datasets are all exported, and their joint fine-mapping runs. */
    private static void collectRuns(Set<String> ids, Map<File, String> out, List<String> jobs) {
        File[] dirs = new File("output/multi_locus").listFiles(File::isDirectory);
        if (dirs != null) for (File d : dirs) {
            File jf = new File(d, "job.json");
            if (!jf.isFile()) continue;
            try {
                Map<String, Object> m = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(jf.toPath()), StandardCharsets.UTF_8)));
                List<Object> pids = m.get("project_ids") instanceof List ? MiniJson.asArray(m.get("project_ids")) : new ArrayList<>();
                if (pids.isEmpty() || !"done".equals(MiniJson.getStr(m, "status", ""))) continue;
                boolean all = true;
                for (Object o : pids) all &= ids.contains(String.valueOf(o));
                if (!all) continue;
                jobs.add(d.getName());
                addTree(d, "output/multi_locus/" + d.getName(), out, f -> true);
            } catch (Exception ignored) {}
        }
        File[] mv = new File("output/mvsusie").listFiles(File::isDirectory);
        if (mv != null && !jobs.isEmpty()) for (File d : mv) {
            File meta = new File(d, "meta.json");
            if (!meta.isFile()) continue;
            try {
                Map<String, Object> m = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(meta.toPath()), StandardCharsets.UTF_8)));
                if (jobs.contains(MiniJson.getStr(m, "job", ""))) addTree(d, "output/mvsusie/" + d.getName(), out, MVSUSIE_KEEP::contains);
            } catch (Exception ignored) {}
        }
    }

    /** The newest run of each tool, chosen as LocalServer.latestRuns does. */
    static List<File> latestRuns(File runsDir) {
        Map<String, Object[]> best = new HashMap<>();
        File[] dirs = runsDir.listFiles(File::isDirectory);
        if (dirs != null) for (File d : dirs) {
            File prov = new File(d, "provenance.json");
            if (!prov.isFile() || !new File(d, "result.tsv").isFile()) continue;
            try {
                Map<String, Object> m = MiniJson.asObject(MiniJson.parse(new String(Files.readAllBytes(prov.toPath()), StandardCharsets.UTF_8)));
                String tool = MiniJson.getStr(m, "tool", "");
                long ts = m.get("timestamp") instanceof Number ? ((Number) m.get("timestamp")).longValue() : d.lastModified();
                Object[] b = best.get(tool);
                if (b == null || ts > (Long) b[0]) best.put(tool, new Object[]{ts, d});
            } catch (Exception ignored) {}
        }
        List<File> out = new ArrayList<>();
        for (Object[] b : best.values()) out.add((File) b[1]);
        return out;
    }

    private static void addTree(File root, String arc, Map<File, String> out, java.util.function.Predicate<String> keep) {
        File[] items = root.listFiles();
        if (items == null) return;
        for (File f : items) {
            if (f.isDirectory()) addTree(f, arc + "/" + f.getName(), out, keep);
            else if (keep.test(f.getName())) out.put(f, arc + "/" + f.getName());
        }
    }

    /** The spellings an absolute path takes in config, JSON and logs: C:/x, C:\x, C:\\x. */
    private static void addRule(List<String[]> rules, String path, String placeholder) {
        String fwd = path.replace('\\', '/');
        while (fwd.endsWith("/")) fwd = fwd.substring(0, fwd.length() - 1);
        if (fwd.length() < 4) return;                                    // never a drive root
        for (String v : new String[]{fwd, fwd.replace("/", "\\"), fwd.replace("/", "\\\\")}) {
            boolean dup = false;
            for (String[] r : rules) dup |= r[0].equals(v);
            if (!dup) rules.add(new String[]{v, placeholder});
        }
    }

    private static boolean isText(String arc) {
        int dot = arc.lastIndexOf('.');
        return arc.endsWith("/config.properties") || (dot > arc.lastIndexOf('/') && TEXT_EXT.contains(arc.substring(dot + 1).toLowerCase(Locale.ROOT)));
    }

    /** A project's project.json: results only unless its summary statistics travel with it. */
    private static String projectOverride(String arc, String text, List<String> ids, boolean includeGwas) {
        if (arc.startsWith("projects/") && arc.endsWith("/config.properties")) {
            // paths in the pack must work on any system: forward slashes after the placeholders
            StringBuilder b = new StringBuilder();
            for (String line : text.split("\n", -1)) {
                int at = line.indexOf("@@LYNX_");
                b.append(at < 0 ? line : line.substring(0, at) + line.substring(at).replace("\\\\", "/").replace('\\', '/')).append('\n');
            }
            text = b.substring(0, b.length() - 1);
        }
        for (String id : ids) {
            if (!arc.equals("projects/" + id + "/project.json")) continue;
            boolean hasGwas = false;
            try {
                Config cfg = Config.loadFromProject(new File("projects", id).getAbsolutePath());
                hasGwas = includeGwas && new File(cfg.gwasFile).isFile();
            } catch (Exception ignored) {}
            ProjectMetadata pm = ProjectMetadata.load(new File("projects", id).getAbsolutePath());
            boolean resultsOnly = !hasGwas || (pm != null && pm.resultsOnly);
            text = text.replaceAll(",?\\s*\"results_only\"\\s*:\\s*(true|false)", "");
            if (resultsOnly) text = text.replaceFirst("\\s*}\\s*$", ",\n  \"results_only\": true\n}\n");
        }
        return text;
    }

    private static void checkCancel(ExampleDataInstaller.Progress p) {
        if (p.cancelled) throw new ExampleDataInstaller.CancelledException();
    }
}
