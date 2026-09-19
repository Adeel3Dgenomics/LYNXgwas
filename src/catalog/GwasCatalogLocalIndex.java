package catalog;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/**
 * Local-first lookup against the public EBI GWAS Catalog bulk associations dump
 * (https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip,
 * ~741MB unzipped, ~1.19M rows). Builds two on-disk indexes over that file in one streaming pass:
 *
 *   1. rsID -> byte offset(s)        — exact-variant lookup (see {@link #lookup}), used for the
 *                                       per-locus "known associations for this SNP" panel.
 *   2. chr -> sorted (pos, offset)[] — region-overlap lookup (see {@link #overlap}), used for the
 *                                       known-vs-novel-*locus* verdict: does ANY previously reported
 *                                       association fall inside this locus's genomic window, whether
 *                                       or not it's the exact same SNP as our own top hit.
 *
 * Both answer with a handful of RandomAccessFile seeks — no network call after the initial download.
 * (The live REST API only gives a p-value inline per association; trait + study each need a separate
 * follow-up GET, which turns a popular-SNP lookup into hundreds of round trips — confirmed hanging
 * past 90s on rs429358/APOE. This file exists specifically to avoid that.)
 */
public class GwasCatalogLocalIndex {

    // Column indices in gwas-catalog-download-associations-alt-full.tsv (0-based)
    private static final int COL_PUBMEDID = 1, COL_LINK = 5, COL_DISEASE_TRAIT = 7, COL_INITIAL_SAMPLE = 8,
        COL_CHR_ID = 11, COL_CHR_POS = 12, COL_MAPPED_GENE = 14, COL_SNPS = 21,
        COL_PVALUE = 27, COL_PVALUE_MLOG = 28, COL_MAPPED_TRAIT = 34, COL_STUDY_ACCESSION = 36;

    private final Path tsvPath;
    private final Path rsidIndexPath;
    private final Path regionIndexPath;
    private Map<String, long[]> rsidIndex;              // rsid -> offsets
    private Map<String, long[]> regionPos;              // chr -> sorted positions
    private Map<String, long[]> regionOff;              // chr -> parallel offsets (same order as regionPos)
    private Map<String, String[]> regionTraitBlob;      // chr -> parallel lowercased trait text (same order)

    public GwasCatalogLocalIndex(String tsvPath, String indexDir) {
        this.tsvPath = Paths.get(tsvPath);
        this.rsidIndexPath = Paths.get(indexDir, "rsid_index.tsv");
        this.regionIndexPath = Paths.get(indexDir, "region_index.tsv");
    }

    public static class Hit {
        public String trait, mappedTrait, pubmedId, studyAccession, initialSampleSize, mappedGene, link;
        public String chr; public long pos;
        public double pvalue;     // 0.0 for anything below ~4.9E-324 — Double can't represent it at all
        public double pvalueMlog; // -log10(p) from the Catalog's own PVALUE_MLOG column; NaN if absent.
                                   // Survives the underflow above since e.g. -log10(6E-363)=362.2 is an
                                   // ordinary double, so this is the only reliable source of magnitude
                                   // once pvalue itself has flattened to zero.

        /** Human-readable p-value, correct even when pvalue underflowed to exactly 0.0. */
        public String pvalueDisplay() {
            if (pvalue > 0 && Double.isFinite(pvalue)) {
                String s = String.format("%.2e", pvalue);
                return s.replaceFirst("e([+-])0*(\\d)", "e$1$2"); // "3.00e-08" -> "3.00e-8"
            }
            if (Double.isFinite(pvalueMlog) && pvalueMlog > 0) {
                return "<1e-" + Math.round(pvalueMlog); // exact mantissa is unrecoverable once underflowed
            }
            return pvalue == 0 ? "0" : "NA";
        }

        @Override public String toString() {
            return String.format("%-55s p=%s  gene=%-12s study=%s  PMID=%s", trait, pvalueDisplay(), mappedGene, studyAccession, pubmedId);
        }
    }

    public enum Verdict { NOVEL, KNOWN_SAME_TRAIT, KNOWN_OTHER_TRAIT, KNOWN_UNSPECIFIED }

    public static class RegionResult {
        public Verdict verdict;
        public List<Hit> hits;   // capped at MAX_OVERLAP_HITS — see classify()
        public int totalCount;   // exact, uncapped — a free by-product of the binary search, no seeks needed
    }

    /** Hard cap on how many hits overlap() actually seeks and parses. A handful of genes (APOE, the
     *  MHC region, etc.) have tens of thousands of reported associations — reading and JSON-serializing
     *  all of them turned one locus's novelty check into a 12-second, 4.5MB response. totalCount above
     *  is still exact (computed from the sorted-position array alone, no seeks), so the known/novel
     *  verdict and the true count are never wrong; only the displayed hit *list* is capped. */
    public static final int MAX_OVERLAP_HITS = 200;

    /** One region-index row: position + byte offset (for capped, on-demand full parsing) + a small
     *  lowercased trait blob kept resident in memory so trait-keyword matching (see classify()) never
     *  needs a disk seek, however many thousands of rows a locus overlaps. */
    private static class RegionEntry {
        long pos, offset;
        String traitBlob;
    }

    /** Builds (or rebuilds) both on-disk indexes. One streaming pass over the full TSV. */
    public void buildIndex() throws IOException {
        System.out.println("[GwasCatalogLocalIndex] Building indexes from " + tsvPath + " ...");
        long t0 = System.currentTimeMillis();
        Map<String, List<Long>> rsidMap = new HashMap<>(2_000_000);
        Map<String, List<RegionEntry>> regionMap = new HashMap<>(); // chr -> entries
        int rows = 0, skippedRows = 0;

        // Raw byte scan, not BufferedReader.readLine() — the file uses CRLF, and byte-counting off a
        // String's length instead of the actual stream undercounts by 1 byte/line and drifts every
        // offset after the first line. Reading raw bytes and finding '\n' ourselves keeps offsets exact.
        try (InputStream in = new BufferedInputStream(new FileInputStream(tsvPath.toFile()), 1 << 20)) {
            ByteArrayOutputStream lineBuf = new ByteArrayOutputStream(4096);
            long offset = 0, lineStart = 0;
            boolean firstLine = true;
            int b;
            while ((b = in.read()) != -1) {
                offset++;
                if (b == '\n') {
                    String line = stripTrailingCr(lineBuf.toString(StandardCharsets.UTF_8));
                    if (!firstLine) {
                        // A single malformed row must never take down a 30s index build — skip and count it.
                        try {
                            String[] cols = line.split("\t", -1);
                            if (cols.length > COL_SNPS) {
                                for (String rsid : cols[COL_SNPS].split("\\s*,\\s*")) {
                                    rsid = rsid.trim();
                                    if (!rsid.isEmpty())
                                        rsidMap.computeIfAbsent(rsid, k -> new ArrayList<>()).add(lineStart);
                                }
                                rows++;
                            }
                            if (cols.length > COL_CHR_POS) {
                                String chr = normalizeChr(firstToken(cols[COL_CHR_ID]));
                                long pos = parseLongOrMinus1(firstToken(cols[COL_CHR_POS]));
                                if (chr != null && pos >= 0) {
                                    RegionEntry re = new RegionEntry();
                                    re.pos = pos;
                                    re.offset = lineStart;
                                    re.traitBlob = traitBlob(cols);
                                    regionMap.computeIfAbsent(chr, k -> new ArrayList<>()).add(re);
                                }
                            }
                        } catch (RuntimeException e) {
                            skippedRows++;
                        }
                    }
                    firstLine = false;
                    lineBuf.reset();
                    lineStart = offset;
                    continue;
                }
                lineBuf.write(b);
            }
        }

        try (BufferedWriter bw = Files.newBufferedWriter(rsidIndexPath, StandardCharsets.UTF_8)) {
            for (Map.Entry<String, List<Long>> e : rsidMap.entrySet()) {
                bw.write(e.getKey()); bw.write('\t');
                StringBuilder sb = new StringBuilder();
                for (long off : e.getValue()) { if (sb.length() > 0) sb.append(','); sb.append(off); }
                bw.write(sb.toString()); bw.newLine();
            }
        }

        try (BufferedWriter bw = Files.newBufferedWriter(regionIndexPath, StandardCharsets.UTF_8)) {
            for (Map.Entry<String, List<RegionEntry>> e : regionMap.entrySet()) {
                List<RegionEntry> entries = e.getValue();
                entries.sort(Comparator.comparingLong(re -> re.pos));
                StringBuilder posSb = new StringBuilder(), offSb = new StringBuilder(), traitSb = new StringBuilder();
                for (RegionEntry re : entries) {
                    if (posSb.length() > 0) { posSb.append(','); offSb.append(','); traitSb.append(TRAIT_SEP); }
                    posSb.append(re.pos); offSb.append(re.offset); traitSb.append(re.traitBlob);
                }
                bw.write(e.getKey()); bw.write('\t');
                bw.write(posSb.toString()); bw.write('\t');
                bw.write(offSb.toString()); bw.write('\t');
                bw.write(traitSb.toString()); bw.newLine();
            }
        }

        System.out.printf("[GwasCatalogLocalIndex] Indexed %d rows (%d skipped), %d distinct rsIDs, %d chromosomes in %.1fs%n",
            rows, skippedRows, rsidMap.size(), regionMap.size(), (System.currentTimeMillis() - t0) / 1000.0);
    }

    // Splits on a real multi-value separator (";", or " x "/" X " for interaction-study notation like
    // "6 x 12") — NOT on a bare "x"/"X" character, since that's also how the Catalog spells the X
    // chromosome itself; splitting on it there produced a string of nothing but delimiters, which
    // String.split() collapses to a zero-length array (trailing empties are stripped) and crashed
    // every request that touched an X-chromosome row.
    /** Control character used to join per-row trait blobs in the region index file. Real trait text
     *  never contains it, unlike comma/pipe/semicolon which do show up in compound trait names. */
    private static final char TRAIT_SEP = '\u0001';

    private static String traitBlob(String[] cols) {
        String disease = cols.length > COL_DISEASE_TRAIT ? cols[COL_DISEASE_TRAIT] : "";
        String mapped = cols.length > COL_MAPPED_TRAIT ? cols[COL_MAPPED_TRAIT] : "";
        String blob = (disease + " " + mapped).toLowerCase(Locale.ROOT).replace(TRAIT_SEP, ' ');
        return blob.trim();
    }

    private static String firstToken(String s) {
        if (s == null || s.isEmpty()) return null;
        String[] parts = s.split("\\s*;\\s*|\\s+[xX]\\s+");
        String t = (parts.length > 0 ? parts[0] : s).trim();
        return t.isEmpty() ? null : t;
    }

    private static long parseLongOrMinus1(String s) {
        if (s == null) return -1;
        try { return Long.parseLong(s.trim()); } catch (NumberFormatException e) { return -1; }
    }

    private static String normalizeChr(String chr) {
        if (chr == null) return null;
        chr = chr.trim();
        if (chr.regionMatches(true, 0, "chr", 0, 3)) chr = chr.substring(3);
        return chr.isEmpty() ? null : chr;
    }

    private static String stripTrailingCr(String s) {
        return s.endsWith("\r") ? s.substring(0, s.length() - 1) : s;
    }

    /** Loads both (already-built) indexes into memory. Call once before lookup()/overlap(). */
    public void loadIndex() throws IOException {
        long t0 = System.currentTimeMillis();

        rsidIndex = new HashMap<>(2_000_000);
        try (BufferedReader br = Files.newBufferedReader(rsidIndexPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                int tab = line.indexOf('\t');
                if (tab < 0) continue;
                rsidIndex.put(line.substring(0, tab), parseLongCsv(line.substring(tab + 1)));
            }
        }

        regionPos = new HashMap<>();
        regionOff = new HashMap<>();
        regionTraitBlob = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(regionIndexPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split("\t", -1);
                if (parts.length != 4) continue;
                regionPos.put(parts[0], parseLongCsv(parts[1]));
                regionOff.put(parts[0], parseLongCsv(parts[2]));
                regionTraitBlob.put(parts[0], parts[3].split(String.valueOf(TRAIT_SEP), -1));
            }
        }

        System.out.printf("[GwasCatalogLocalIndex] Loaded %d rsIDs + %d chromosomes' regions into memory in %.1fs%n",
            rsidIndex.size(), regionPos.size(), (System.currentTimeMillis() - t0) / 1000.0);
    }

    private static long[] parseLongCsv(String s) {
        if (s.isEmpty()) return new long[0];
        String[] parts = s.split(",");
        long[] out = new long[parts.length];
        for (int i = 0; i < parts.length; i++) out[i] = Long.parseLong(parts[i]);
        return out;
    }

    public boolean isLoaded() { return rsidIndex != null; }

    /** Instant exact-rsID lookup — no network. Empty list, not null, if the rsID isn't in the snapshot. */
    public List<Hit> lookup(String rsid) throws IOException {
        if (rsidIndex == null) throw new IllegalStateException("call loadIndex() first");
        long[] offsets = rsidIndex.get(rsid);
        if (offsets == null) return new ArrayList<>();
        try (RandomAccessFile raf = new RandomAccessFile(tsvPath.toFile(), "r")) {
            List<Hit> hits = new ArrayList<>();
            for (long off : offsets) hits.add(readHitAt(raf, off));
            hits.sort(Comparator.comparingDouble(h -> h.pvalue));
            return hits;
        }
    }

    /**
     * Every previously reported association whose position falls within [start, end] on this
     * chromosome — regardless of whether it's the same rsID as our own top SNP. This is the check
     * that actually answers "is this locus known", since a locus's originally-reported tag SNP can
     * differ from ours by a few kb and still represent the same signal.
     */
    public List<Hit> overlap(String chr, long start, long end) throws IOException {
        return overlap(chr, start, end, MAX_OVERLAP_HITS);
    }

    /** Like overlap(), but seeks+parses at most maxHits rows (in position order) regardless of how
     *  many actually fall in the region — see MAX_OVERLAP_HITS. Use overlapCount() for the exact total. */
    public List<Hit> overlap(String chr, long start, long end, int maxHits) throws IOException {
        if (regionPos == null) throw new IllegalStateException("call loadIndex() first");
        chr = normalizeChr(chr);
        long[] positions = regionPos.get(chr);
        long[] offsets = regionOff.get(chr);
        List<Hit> hits = new ArrayList<>();
        if (positions == null || positions.length == 0) return hits;

        int lo = lowerBound(positions, start);
        if (lo >= positions.length || positions[lo] > end) return hits;

        int limit = Math.min(positions.length, lo + maxHits);
        try (RandomAccessFile raf = new RandomAccessFile(tsvPath.toFile(), "r")) {
            for (int i = lo; i < limit && positions[i] <= end; i++) {
                hits.add(readHitAt(raf, offsets[i]));
            }
        }
        hits.sort(Comparator.comparingDouble(h -> h.pvalue));
        return hits;
    }

    /** Exact count of associations overlapping this region — pure binary search, no disk seeks at all,
     *  so this stays cheap even for a locus with tens of thousands of hits. */
    public int overlapCount(String chr, long start, long end) {
        if (regionPos == null) throw new IllegalStateException("call loadIndex() first");
        chr = normalizeChr(chr);
        long[] positions = regionPos.get(chr);
        if (positions == null || positions.length == 0) return 0;
        int lo = lowerBound(positions, start);
        int hi = lo;
        while (hi < positions.length && positions[hi] <= end) hi++;
        return hi - lo;
    }

    /** Convenience: overlap() + a novel/known-same-trait/known-other-trait verdict. traitKeyword may be null. */
    public RegionResult classify(String chr, long start, long end, String traitKeyword) throws IOException {
        RegionResult r = new RegionResult();
        r.totalCount = overlapCount(chr, start, end);
        r.hits = r.totalCount == 0 ? new ArrayList<>() : overlap(chr, start, end);
        if (r.hits.isEmpty()) {
            r.verdict = Verdict.NOVEL;
        } else if (traitKeyword == null || traitKeyword.trim().isEmpty()) {
            r.verdict = Verdict.KNOWN_UNSPECIFIED;
        } else {
            // Deliberately NOT r.hits (capped at MAX_OVERLAP_HITS, sorted by significance) — a locus
            // like APOE has 13,000+ associations, and the 2 that happen to mention a rarer trait can
            // easily rank outside the top 200 by p-value while still being perfectly real matches.
            // regionTraitBlob covers every overlapping row with zero disk seeks, so the trait-match
            // verdict stays exact regardless of how large the capped display list is.
            boolean sameTrait = containsTraitKeyword(chr, start, end, traitKeyword.toLowerCase(Locale.ROOT));
            r.verdict = sameTrait ? Verdict.KNOWN_SAME_TRAIT : Verdict.KNOWN_OTHER_TRAIT;
        }
        return r;
    }

    /** True if any association overlapping this region has traitKeywordLower as a substring of its
     *  disease-trait or mapped-trait text. Pure in-memory scan (regionTraitBlob), no disk I/O — exact
     *  even for a locus with tens of thousands of overlapping rows. */
    private boolean containsTraitKeyword(String chr, long start, long end, String traitKeywordLower) {
        chr = normalizeChr(chr);
        long[] positions = regionPos.get(chr);
        String[] traits = regionTraitBlob.get(chr);
        if (positions == null || traits == null) return false;
        int lo = lowerBound(positions, start);
        for (int i = lo; i < positions.length && positions[i] <= end; i++) {
            if (traits[i].contains(traitKeywordLower)) return true;
        }
        return false;
    }

    /** First index i such that positions[i] >= target (binary search over a sorted ascending array). */
    private static int lowerBound(long[] positions, long target) {
        int lo = 0, hi = positions.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (positions[mid] < target) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    private Hit readHitAt(RandomAccessFile raf, long offset) throws IOException {
        raf.seek(offset);
        String line = readUtf8Line(raf);
        Hit h = new Hit();
        if (line == null) return h;
        String[] c = line.split("\t", -1);
        if (c.length <= COL_STUDY_ACCESSION) return h;
        h.trait = c[COL_DISEASE_TRAIT];
        h.mappedTrait = c[COL_MAPPED_TRAIT];
        h.pubmedId = c[COL_PUBMEDID];
        h.studyAccession = c[COL_STUDY_ACCESSION];
        h.initialSampleSize = c[COL_INITIAL_SAMPLE];
        h.mappedGene = c[COL_MAPPED_GENE];
        h.link = c[COL_LINK];
        h.chr = firstToken(c[COL_CHR_ID]);
        h.pos = parseLongOrMinus1(firstToken(c[COL_CHR_POS]));
        try { h.pvalue = Double.parseDouble(c[COL_PVALUE]); } catch (NumberFormatException e) { h.pvalue = Double.NaN; }
        h.pvalueMlog = Double.NaN;
        if (c.length > COL_PVALUE_MLOG) {
            try { h.pvalueMlog = Double.parseDouble(c[COL_PVALUE_MLOG]); } catch (NumberFormatException ignored) {}
        }
        return h;
    }

    private static String readUtf8Line(RandomAccessFile raf) throws IOException {
        ByteArrayOutputStream buf = new ByteArrayOutputStream(256);
        int b;
        while ((b = raf.read()) != -1 && b != '\n') buf.write(b);
        if (b == -1 && buf.size() == 0) return null;
        return stripTrailingCr(buf.toString(StandardCharsets.UTF_8));
    }
}
