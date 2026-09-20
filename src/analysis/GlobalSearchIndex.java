import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.time.*;
import java.time.format.DateTimeFormatter;
import java.util.*;

/**
 * Process-lifetime, lazily-built, in-memory index over every project under {@code projects/*}
 * (any directory with a {@code data/manifest.json}), answering "which of my projects have a
 * significant hit at gene X / SNP rsY / position chr:pos" without opening each project one at a
 * time.
 *
 * Only SNPs reaching a significance threshold (default {@link #DEFAULT_THRESHOLD}, 5e-8) are ever
 * indexed — the threshold is a build-time parameter (like the rest of this app's simple,
 * no-file-watcher, build-once-and-cache philosophy — see {@link catalog.GwasCatalogLocalIndex} for
 * the same shape of index elsewhere in this codebase, though that one persists to disk and this one
 * deliberately doesn't). A search at a different threshold than the currently cached index rebuilds
 * it once for that threshold (see {@link #getOrBuild}) and reuses it after that.
 *
 * Gene matches are always tagged honestly as "contains" (the SNP's position falls inside that
 * gene's own [start,end] in that locus's "genes" array) or "nearest" (no gene in the locus contains
 * it, so it's recorded against every name in the locus's manifest-level nearest_genes instead) —
 * the two are never blurred together.
 *
 * Position search is build-aware:
 *   - query build == project build            -> direct comparison, no liftover.
 *   - query GRCh38, project GRCh37            -> project's SNP positions are lifted once via
 *                                                 {@link GenomeLiftover#toGRCh38} at build time and
 *                                                 cached (never re-lifted per query).
 *   - query GRCh37, project GRCh38            -> NOT SUPPORTED (no reverse liftover exists in this
 *                                                 codebase). Such projects are never silently
 *                                                 compared or omitted — they come back in the
 *                                                 response's notCompared list instead.
 */
public class GlobalSearchIndex {

    public static final double DEFAULT_THRESHOLD = 5e-8;

    // ── Process-lifetime singleton, used by LocalServer ─────────────────────
    private static volatile GlobalSearchIndex INSTANCE;
    private static volatile double instanceThreshold = Double.NaN;

    /** Builds the index on first call (or if the requested threshold differs from the cached
     *  build's), otherwise returns the cached instance. */
    public static synchronized GlobalSearchIndex getOrBuild(File projectsRoot, double threshold) {
        if (INSTANCE == null || Double.compare(instanceThreshold, threshold) != 0) {
            INSTANCE = build(projectsRoot, threshold);
            instanceThreshold = threshold;
        }
        return INSTANCE;
    }

    /** Forces a fresh build regardless of what's cached. Powers POST /api/search/rebuild-index. */
    public static synchronized GlobalSearchIndex forceRebuild(File projectsRoot, double threshold) {
        INSTANCE = build(projectsRoot, threshold);
        instanceThreshold = threshold;
        return INSTANCE;
    }

    // ── Summary stats from the last build (logged to stdout too) ───────────
    public int projectCount, lociCount, snpCount;
    public long buildMillis;
    public String builtAt = "";
    public double threshold = DEFAULT_THRESHOLD;

    // ── Indexes ──────────────────────────────────────────────────────────
    private final Map<String, List<Hit>> byGeneUpper = new HashMap<>();
    private final Map<String, List<Hit>> byRsid       = new HashMap<>();
    // Native (un-lifted) position index, keyed by the project's own declared genome.build, then
    // by "chr:pos". Kept split per build so two projects on different builds can never collide on
    // a coincidentally-identical raw chr:pos pair.
    private final Map<String, Map<String, List<Hit>>> byPosNative = new HashMap<>();
    // GRCh38-equivalent position index: native GRCh38 hits, plus GRCh37 hits successfully lifted
    // (cached at build time, never re-lifted per query).
    private final Map<String, List<Hit>> byPosGRCh38 = new HashMap<>();
    private final List<ProjectBuildInfo> projectBuilds = new ArrayList<>();

    public static class ProjectBuildInfo {
        public String id, name, build;
    }

    /** One result row. gene_name/match_type are null for rsID/position search results — those
     *  searches don't compute gene overlap, so they honestly report "unknown" rather than guessing. */
    public static class Hit {
        public String projectId = "", projectName = "";
        public int    locusIndex;
        public String locusName = "";
        public String geneName;      // nullable
        public String matchType;     // "contains" | "nearest" | null
        public String snpId = "";
        public String chr = "";
        public long   pos;
        public double pvalue;
        public Double effectValue;   // nullable
        public String effectType = ""; // "beta" | "OR" | "logOR" | "" (unknown)

        public String toJson() {
            StringBuilder j = new StringBuilder();
            j.append('{');
            kv(j, "project_id", projectId, true);
            kv(j, "project_name", projectName, false);
            j.append(",\"locus_index\":").append(locusIndex);
            kv(j, "locus_name", locusName, false);
            j.append(",\"gene_name\":").append(geneName == null ? "null" : quoted(geneName));
            j.append(",\"match_type\":").append(matchType == null ? "null" : quoted(matchType));
            kv(j, "snp_id", snpId, false);
            kv(j, "chr", chr, false);
            j.append(",\"pos\":").append(pos);
            j.append(",\"pvalue\":").append(num(pvalue));
            j.append(",\"effect_value\":").append(effectValue == null ? "null" : num(effectValue));
            kv(j, "effect_type", effectType, false);
            j.append('}');
            return j.toString();
        }
    }

    /** A project that could not be honestly compared for a given position query (build mismatch
     *  with no supported liftover direction) — surfaced explicitly, never silently dropped. */
    public static class NotCompared {
        public String projectId, projectName, reason;

        public String toJson() {
            StringBuilder j = new StringBuilder();
            j.append('{');
            kv(j, "project_id", projectId, true);
            kv(j, "project_name", projectName, false);
            kv(j, "reason", reason, false);
            j.append('}');
            return j.toString();
        }
    }

    public static class Result {
        public List<Hit> hits = new ArrayList<>();
        public List<NotCompared> notCompared = new ArrayList<>();
    }

    // ── Search ───────────────────────────────────────────────────────────

    public Result searchGene(String geneName) {
        Result r = new Result();
        if (geneName == null || geneName.trim().isEmpty()) return r;
        List<Hit> hits = byGeneUpper.get(geneName.trim().toUpperCase(Locale.ROOT));
        if (hits != null) r.hits.addAll(hits);
        return r;
    }

    public Result searchSnp(String rsid) {
        Result r = new Result();
        if (rsid == null || rsid.trim().isEmpty()) return r;
        List<Hit> hits = byRsid.get(rsid.trim());
        if (hits != null) r.hits.addAll(hits);
        return r;
    }

    /** queryBuild: "GRCh37" or "GRCh38" (defaults to "GRCh37" if blank, matching Config's own default). */
    public Result searchPosition(String chr, long pos, String queryBuild) {
        Result r = new Result();
        String key = normalizeChr(chr) + ":" + pos;
        String qb = (queryBuild == null || queryBuild.trim().isEmpty()) ? "GRCh37" : queryBuild.trim();

        if ("GRCh38".equalsIgnoreCase(qb)) {
            // Native GRCh38 projects + GRCh37 projects successfully lifted are already merged here.
            List<Hit> hits = byPosGRCh38.get(key);
            if (hits != null) r.hits.addAll(hits);
        } else if ("GRCh37".equalsIgnoreCase(qb)) {
            Map<String, List<Hit>> nativeGrch37 = byPosNative.get("GRCh37");
            List<Hit> hits = nativeGrch37 != null ? nativeGrch37.get(key) : null;
            if (hits != null) r.hits.addAll(hits);
            // No GRCh38->GRCh37 liftover exists anywhere in this codebase (GenomeLiftover is
            // GRCh37->GRCh38 only). Every GRCh38 project is therefore explicitly flagged as not
            // compared for a GRCh37 query, rather than silently skipped or wrongly compared.
            for (ProjectBuildInfo pbi : projectBuilds) {
                if ("GRCh38".equalsIgnoreCase(pbi.build)) {
                    NotCompared nc = new NotCompared();
                    nc.projectId = pbi.id;
                    nc.projectName = pbi.name;
                    nc.reason = "not compared (GRCh38→GRCh37 liftover unavailable)";
                    r.notCompared.add(nc);
                }
            }
        } else {
            // Unrecognized build string on the query itself — fall back to a same-build-only
            // comparison rather than guessing a liftover direction.
            Map<String, List<Hit>> native_ = byPosNative.get(qb);
            List<Hit> hits = native_ != null ? native_.get(key) : null;
            if (hits != null) r.hits.addAll(hits);
        }
        sortHits(r.hits);
        return r;
    }

    public String summaryJson() {
        StringBuilder j = new StringBuilder();
        j.append('{');
        j.append("\"project_count\":").append(projectCount).append(',');
        j.append("\"loci_count\":").append(lociCount).append(',');
        j.append("\"snp_count\":").append(snpCount).append(',');
        j.append("\"build_ms\":").append(buildMillis).append(',');
        j.append("\"threshold\":").append(num(threshold)).append(',');
        kv(j, "built_at", builtAt, true); // "first" here just means "no comma" — one was already appended above
        j.append('}');
        return j.toString();
    }

    // ── Build ────────────────────────────────────────────────────────────

    private static class SnpRecord {
        String projectId, projectName, locusName, snpId, chr, effectType;
        int locusIndex;
        long pos;
        double pvalue;
        Double effectValue;
    }

    public static GlobalSearchIndex build(File projectsRoot, double threshold) {
        long t0 = System.currentTimeMillis();
        GlobalSearchIndex idx = new GlobalSearchIndex();
        idx.threshold = threshold;

        int projects = 0, lociTotal = 0, snpsTotal = 0;

        File[] children = projectsRoot.isDirectory() ? projectsRoot.listFiles() : null;
        if (children != null) {
            Arrays.sort(children, Comparator.comparing(File::getName));
            for (File dir : children) {
                if (!dir.isDirectory()) continue;
                File manifestFile = new File(dir, "data/manifest.json");
                if (!manifestFile.isFile()) continue;

                String projectId = dir.getName();
                try {
                    Config cfg;
                    try {
                        cfg = Config.loadFromProject(dir.getAbsolutePath());
                    } catch (Exception e) {
                        System.err.println("[GlobalSearchIndex] Skipping '" + projectId + "' (config): " + e.getMessage());
                        continue;
                    }
                    projects++;

                    ProjectMetadata pm = ProjectMetadata.load(dir.getAbsolutePath());
                    String projectName = (pm != null && pm.name != null && !pm.name.isEmpty()) ? pm.name : projectId;
                    String build = (cfg.genomeBuild == null || cfg.genomeBuild.trim().isEmpty()) ? "GRCh37" : cfg.genomeBuild.trim();
                    String effectType = cfg.effectType == null ? "" : cfg.effectType.trim();

                    ProjectBuildInfo pbi = new ProjectBuildInfo();
                    pbi.id = projectId; pbi.name = projectName; pbi.build = build;
                    idx.projectBuilds.add(pbi);

                    String manifestJson = readFile(manifestFile);
                    List<Integer> lociIndices = new ArrayList<>();
                    for (String locusObj : extractObjectArray(manifestJson, "loci")) {
                        long li = extractLong(locusObj, "index", -1);
                        if (li >= 0) lociIndices.add((int) li);
                    }

                    for (int locusIdx : lociIndices) {
                        File locusFile = new File(dir, "data/locus_" + locusIdx + ".json");
                        if (!locusFile.isFile()) continue;
                        String lj;
                        try { lj = readFile(locusFile); }
                        catch (IOException e) { continue; }
                        lociTotal++;

                        String locusName = extractStr(lj, "locus_name");
                        if (locusName == null) locusName = "Locus " + locusIdx;
                        String locusChr = extractStr(lj, "chr");

                        List<String> nearestGenes = extractStringArray(lj, "nearest_genes");

                        List<String> geneNames = new ArrayList<>();
                        List<long[]> geneRanges = new ArrayList<>();
                        for (String geneObj : extractObjectArray(lj, "genes")) {
                            String gn = extractStr(geneObj, "gene_name");
                            if (gn == null) continue;
                            long gs = extractLong(geneObj, "start", Long.MIN_VALUE);
                            long ge = extractLong(geneObj, "end", Long.MAX_VALUE);
                            geneNames.add(gn);
                            geneRanges.add(new long[]{gs, ge});
                        }

                        for (String snpObj : extractObjectArray(lj, "gwas_snps")) {
                            double pvalue = extractDouble(snpObj, "pvalue", Double.NaN);
                            if (Double.isNaN(pvalue) || pvalue > threshold) continue; // only significant SNPs are indexed

                            SnpRecord rec = new SnpRecord();
                            rec.projectId = projectId;
                            rec.projectName = projectName;
                            rec.locusIndex = locusIdx;
                            rec.locusName = locusName;
                            String sid = extractStr(snpObj, "id");
                            rec.snpId = sid == null ? "" : sid;
                            String schr = extractStr(snpObj, "chr");
                            rec.chr = schr != null ? schr : (locusChr == null ? "" : locusChr);
                            rec.pos = extractLong(snpObj, "pos", -1);
                            if (rec.pos < 0) continue;
                            rec.pvalue = pvalue;
                            double betaRaw = extractDouble(snpObj, "beta", Double.NaN);
                            rec.effectValue = Double.isNaN(betaRaw) ? null : betaRaw;
                            rec.effectType = effectType;

                            snpsTotal++;

                            // ---- gene index: "contains" if inside any of this locus's own genes,
                            //      else "nearest" fallback against the locus's manifest-level nearest_genes ----
                            boolean contained = false;
                            for (int gi = 0; gi < geneNames.size(); gi++) {
                                long[] range = geneRanges.get(gi);
                                if (rec.pos >= range[0] && rec.pos <= range[1]) {
                                    contained = true;
                                    idx.addGeneHit(geneNames.get(gi), rec, "contains");
                                }
                            }
                            if (!contained) {
                                for (String ng : nearestGenes) idx.addGeneHit(ng, rec, "nearest");
                            }

                            // ---- rsID index ----
                            if (!rec.snpId.isEmpty()) {
                                idx.byRsid.computeIfAbsent(rec.snpId, k -> new ArrayList<>()).add(idx.toHit(rec, null, null));
                            }

                            // ---- position index: native build ----
                            String posKey = normalizeChr(rec.chr) + ":" + rec.pos;
                            idx.byPosNative.computeIfAbsent(build, k -> new HashMap<>())
                                .computeIfAbsent(posKey, k -> new ArrayList<>()).add(idx.toHit(rec, null, null));

                            // ---- position index: GRCh38-effective ----
                            if ("GRCh38".equalsIgnoreCase(build)) {
                                idx.byPosGRCh38.computeIfAbsent(posKey, k -> new ArrayList<>()).add(idx.toHit(rec, null, null));
                            } else if ("GRCh37".equalsIgnoreCase(build)) {
                                GenomeLiftover.Result lifted = GenomeLiftover.toGRCh38(rec.chr, rec.pos, rec.pos);
                                if (lifted.ok) {
                                    String liftedKey = normalizeChr(lifted.chr) + ":" + lifted.start;
                                    idx.byPosGRCh38.computeIfAbsent(liftedKey, k -> new ArrayList<>()).add(idx.toHit(rec, null, null));
                                }
                                // else: chain file missing, or no confident mapping for this exact
                                // position — GenomeLiftover already logs the reason once; this SNP is
                                // simply excluded from the GRCh38-effective index (never silently
                                // compared with un-lifted numbers).
                            }
                        }
                    }
                } catch (Exception e) {
                    System.err.println("[GlobalSearchIndex] Error indexing project '" + projectId + "': " + e);
                }
            }
        }

        for (List<Hit> l : idx.byGeneUpper.values()) sortHits(l);
        for (List<Hit> l : idx.byRsid.values()) sortHits(l);
        for (Map<String, List<Hit>> m : idx.byPosNative.values()) for (List<Hit> l : m.values()) sortHits(l);
        for (List<Hit> l : idx.byPosGRCh38.values()) sortHits(l);

        idx.projectCount = projects;
        idx.lociCount = lociTotal;
        idx.snpCount = snpsTotal;
        idx.buildMillis = System.currentTimeMillis() - t0;
        idx.builtAt = Instant.now().atOffset(ZoneOffset.UTC).format(DateTimeFormatter.ISO_OFFSET_DATE_TIME);

        System.out.printf("[GlobalSearchIndex] Built in %dms: %d project(s), %d loci, %d significant SNP(s) (threshold=%.3g)%n",
            idx.buildMillis, projects, lociTotal, snpsTotal, threshold);

        return idx;
    }

    private void addGeneHit(String geneName, SnpRecord rec, String matchType) {
        if (geneName == null || geneName.trim().isEmpty()) return;
        Hit h = toHit(rec, geneName, matchType);
        byGeneUpper.computeIfAbsent(geneName.trim().toUpperCase(Locale.ROOT), k -> new ArrayList<>()).add(h);
    }

    private Hit toHit(SnpRecord rec, String geneName, String matchType) {
        Hit h = new Hit();
        h.projectId = rec.projectId;
        h.projectName = rec.projectName;
        h.locusIndex = rec.locusIndex;
        h.locusName = rec.locusName;
        h.geneName = geneName;
        h.matchType = matchType;
        h.snpId = rec.snpId;
        h.chr = rec.chr;
        h.pos = rec.pos;
        h.pvalue = rec.pvalue;
        h.effectValue = rec.effectValue;
        h.effectType = rec.effectType;
        return h;
    }

    private static void sortHits(List<Hit> hits) {
        hits.sort(Comparator.comparingDouble(h -> h.pvalue));
    }

    private static String normalizeChr(String chr) {
        if (chr == null) return "";
        chr = chr.trim();
        if (chr.regionMatches(true, 0, "chr", 0, 3)) chr = chr.substring(3);
        return chr;
    }

    private static String readFile(File f) throws IOException {
        return new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
    }

    // ── Minimal hand-rolled JSON helpers (no external dependency — same convention as the rest
    //    of this codebase, e.g. LocalServer's extractStr/extractObjectArray). ──────────────────

    private static String rawValue(String json, String key) {
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

    private static String extractStr(String json, String key) { return rawValue(json, key); }

    private static double extractDouble(String json, String key, double dflt) {
        String v = rawValue(json, key);
        if (v == null || v.isEmpty()) return dflt;
        try { return Double.parseDouble(v); } catch (NumberFormatException e) { return dflt; }
    }

    private static long extractLong(String json, String key, long dflt) {
        String v = rawValue(json, key);
        if (v == null || v.isEmpty()) return dflt;
        try { return Long.parseLong(v); } catch (NumberFormatException e) { return dflt; }
    }

    private static List<String> extractStringArray(String json, String key) {
        List<String> result = new ArrayList<>();
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return result;
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '[') i++;
        if (i >= json.length()) return result;
        int end = json.indexOf(']', i);
        if (end < 0) return result;
        String arrContent = json.substring(i + 1, end);
        for (String part : arrContent.split(",")) {
            part = part.trim();
            if (part.startsWith("\"") && part.endsWith("\""))
                result.add(part.substring(1, part.length() - 1));
        }
        return result;
    }

    /** Splits a top-level JSON array (found via its key) into its object-substrings, one per
     *  top-level element. Works for flat objects (no nested objects), which is all the schemas
     *  this class reads (gwas_snps[], genes[], manifest loci[]) ever contain. */
    private static List<String> extractObjectArray(String json, String key) {
        List<String> result = new ArrayList<>();
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return result;
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '[') i++;
        if (i >= json.length()) return result;
        int arrEnd = i, depth = 0;
        for (int j = i; j < json.length(); j++) {
            char c = json.charAt(j);
            if (c == '[') depth++;
            else if (c == ']') { depth--; if (depth == 0) { arrEnd = j; break; } }
        }
        String inner = json.substring(i + 1, arrEnd);
        int objDepth = 0, objStart = -1;
        for (int j = 0; j < inner.length(); j++) {
            char c = inner.charAt(j);
            if (c == '{') { if (objDepth == 0) objStart = j; objDepth++; }
            else if (c == '}') { objDepth--; if (objDepth == 0 && objStart >= 0) result.add(inner.substring(objStart, j + 1)); }
        }
        return result;
    }

    private static void kv(StringBuilder j, String key, String value, boolean first) {
        if (!first) j.append(',');
        j.append('"').append(key).append("\":\"").append(esc(value)).append('"');
    }

    private static String quoted(String s) { return "\"" + esc(s) + "\""; }

    private static String num(double d) {
        return Double.isNaN(d) || Double.isInfinite(d) ? "null" : String.valueOf(d);
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
