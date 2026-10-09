import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;

/**
 * Maintenance: rewrite the genes and nearest genes of every stored locus with the current gene lookup.
 *
 * Before the GffParser.overlapping fix, a long gene covering a locus was missed whenever a shorter gene sat in
 * between, so stored locus files lack genes (gene track) and some nearest genes are wrong. This recomputes both for
 * each data/locus_N.json (and its .js copy, and the project's manifest) exactly as the pipeline does — genes over
 * the padded window, nearest genes to the locus midpoint — and replaces only those two arrays, leaving every other
 * byte as written. Originals go to <backup>/<project>/ before a file is changed.
 *
 *   java -cp "bin;lib/*" RepairLocusGenes <projects_dir> <backup_dir> [--dry-run]
 */
public class RepairLocusGenes {

    public static void main(String[] args) throws Exception {
        File projects = new File(args[0]);
        File backup = new File(args[1]);
        boolean dry = Arrays.asList(args).contains("--dry-run");
        Map<String, GffParser> gffByPath = new HashMap<>();
        int nLoci = 0, nGenes = 0, nNearest = 0;
        File[] dirs = projects.listFiles(File::isDirectory);
        if (dirs == null) throw new IOException("no projects in " + projects);
        Arrays.sort(dirs);
        for (File pdir : dirs) {
            File data = new File(pdir, "data");
            if (!data.isDirectory() || !new File(pdir, "config.properties").isFile()) continue;
            Config cfg = Config.loadFromProject(pdir.getAbsolutePath());
            if (cfg.gff3File == null || !new File(cfg.gff3File).isFile()) {
                System.out.println("  " + pdir.getName() + ": no gene annotation configured, skipped");
                continue;
            }
            String gp = new File(cfg.gff3File).getAbsolutePath();
            GffParser gff = gffByPath.get(gp);
            if (gff == null) { gff = GffParser.parse(cfg); gffByPath.put(gp, gff); }

            Map<Integer, List<String>> nearestByIndex = new HashMap<>();
            int pLoci = 0, pGenes = 0, pNearest = 0;
            File[] files = data.listFiles((d, n) -> n.matches("locus_\\d+\\.json"));
            if (files == null) continue;
            for (File jf : files) {
                int idx = Integer.parseInt(jf.getName().replaceAll("\\D", ""));
                String txt = new String(Files.readAllBytes(jf.toPath()), StandardCharsets.UTF_8);
                String chr = str(txt, "chr");
                long start = num(txt, "start"), end = num(txt, "end"), ps = num(txt, "padded_start"), pe = num(txt, "padded_end");
                if (chr == null || ps < 0 || pe < 0) continue;
                List<Gene> genes = gff.overlapping(chr.replaceFirst("^chr", ""), ps, pe);
                long mid = (start + end) / 2;
                List<String> nearest = gff.nearestGeneNames(mid, new ArrayList<>(genes));
                String genesJson = JsonExporter.genesToJson(genes);
                String nearestJson = JsonExporter.stringsToJson(nearest);
                String oldGenes = arrayAt(txt, "\"genes\"");
                String oldNearest = arrayAt(txt, "\"nearest_genes\"");
                boolean gDiff = oldGenes != null && !oldGenes.equals(genesJson);
                boolean nDiff = oldNearest != null && !oldNearest.equals(nearestJson);
                nearestByIndex.put(idx, nearest);
                if (!gDiff && !nDiff) continue;
                pLoci++;
                if (gDiff) pGenes += Math.max(0, genes.size() - countGenes(oldGenes));
                if (nDiff) pNearest++;
                if (dry) continue;
                for (File target : new File[]{jf, new File(data, "locus_" + idx + ".js")}) {
                    if (!target.isFile()) continue;
                    String t = new String(Files.readAllBytes(target.toPath()), StandardCharsets.UTF_8);
                    backupOnce(target, new File(backup, pdir.getName()));
                    if (gDiff) t = replaceArray(t, "\"genes\"", genesJson);
                    if (nDiff) t = replaceArray(t, "\"nearest_genes\"", nearestJson);
                    Files.write(target.toPath(), t.getBytes(StandardCharsets.UTF_8));
                }
            }
            // the loci list (manifest.json / manifest.js) carries the nearest genes too
            if (!dry && pNearest > 0) {
                for (File mf : new File[]{new File(data, "manifest.json"), new File(data, "manifest.js")}) {
                    if (!mf.isFile()) continue;
                    String t = new String(Files.readAllBytes(mf.toPath()), StandardCharsets.UTF_8);
                    String u = updateManifest(t, nearestByIndex);
                    if (!u.equals(t)) { backupOnce(mf, new File(backup, pdir.getName())); Files.write(mf.toPath(), u.getBytes(StandardCharsets.UTF_8)); }
                }
            }
            if (pLoci > 0)
                System.out.printf("  %-28s %4d loci changed: +%d genes, %d with new nearest genes%n", pdir.getName(), pLoci, pGenes, pNearest);
            nLoci += pLoci; nGenes += pGenes; nNearest += pNearest;
        }
        System.out.printf("%s: %d loci, +%d genes, %d loci with new nearest genes%n", dry ? "would change" : "changed", nLoci, nGenes, nNearest);
    }

    // ── small JSON helpers working on the text, so nothing else in a file is re-written ──

    private static String str(String t, String key) {
        Matcher m = Pattern.compile("\"" + key + "\":\"([^\"]*)\"").matcher(t);
        return m.find() ? m.group(1) : null;
    }

    private static long num(String t, String key) {
        Matcher m = Pattern.compile("\"" + key + "\":(-?\\d+)").matcher(t);
        return m.find() ? Long.parseLong(m.group(1)) : -1;
    }

    /** [start, end) of the array value that follows the first occurrence of key, or null. */
    private static int[] arraySpan(String t, String key) {
        int k = t.indexOf(key + ":");
        if (k < 0) return null;
        int a = k + key.length() + 1;
        while (a < t.length() && t.charAt(a) == ' ') a++;
        if (a >= t.length() || t.charAt(a) != '[') return null;
        int depth = 0;
        boolean inStr = false;
        for (int i = a; i < t.length(); i++) {
            char c = t.charAt(i);
            if (inStr) { if (c == '\\') i++; else if (c == '"') inStr = false; continue; }
            if (c == '"') inStr = true;
            else if (c == '[') depth++;
            else if (c == ']' && --depth == 0) return new int[]{a, i + 1};
        }
        return null;
    }

    private static String arrayAt(String t, String key) {
        int[] s = arraySpan(t, key);
        return s == null ? null : t.substring(s[0], s[1]);
    }

    private static String replaceArray(String t, String key, String value) {
        int[] s = arraySpan(t, key);
        return s == null ? t : t.substring(0, s[0]) + value + t.substring(s[1]);
    }

    private static int countGenes(String genesArray) {
        Matcher m = Pattern.compile("\"gene_id\":").matcher(genesArray);
        int n = 0;
        while (m.find()) n++;
        return n;
    }

    /** In the manifest, each locus object has "index":N ... "nearest_genes":[...]: replace those arrays. */
    private static String updateManifest(String t, Map<Integer, List<String>> nearest) {
        StringBuilder out = new StringBuilder();
        Matcher m = Pattern.compile("\"index\":(\\d+)").matcher(t);
        int pos = 0;
        while (m.find(pos)) {
            int idx = Integer.parseInt(m.group(1));
            int nextObj = t.indexOf("\"index\":", m.end());
            String seg = t.substring(m.start(), nextObj < 0 ? t.length() : nextObj);
            List<String> ng = nearest.get(idx);
            if (ng != null) {
                int[] s = arraySpan(seg, "\"nearest_genes\"");
                if (s != null) seg = seg.substring(0, s[0]) + JsonExporter.stringsToJson(ng) + seg.substring(s[1]);
            }
            out.append(t, pos, m.start()).append(seg);
            pos = nextObj < 0 ? t.length() : nextObj;
            if (nextObj < 0) break;
        }
        out.append(t.substring(Math.min(pos, t.length())));
        return out.toString();
    }

    private static void backupOnce(File f, File dir) throws IOException {
        dir.mkdirs();
        File b = new File(dir, f.getName());
        if (!b.exists()) Files.copy(f.toPath(), b.toPath());
    }
}
