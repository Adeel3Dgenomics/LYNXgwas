import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Regression test for RegulatoryPeakIndex, using a small hand-written synthetic narrowPeak fixture
 * (real narrowPeak format, made-up positions/signal values chosen so the expected overlap result
 * for each query can be worked out by hand) rather than the real ~200k-line Roadmap files, so this
 * test is fast, deterministic, and independent of whatever real data happens to be downloaded.
 *
 * Fixture for (EID="TESTE", mark="H3K27ac"), query window [1000, 2000] on chr1:
 *   P_before     chr1  500   900   -> end(900)  <  1000 (window start)      -> NOT included
 *   P_edgeStart  chr1  800   1000  -> end(1000) == 1000 (window start)      -> included (boundary)
 *   P_inside     chr1  1200  1300  -> fully inside                          -> included
 *   P_edgeEnd    chr1  2000  2100  -> start(2000) == 2000 (window end)      -> included (boundary)
 *   P_after      chr1  2200  2300  -> start(2200) >  2000 (window end)      -> NOT included
 *   P_otherChr   chr2  1200  1300  -> different chromosome                  -> NOT included
 * Expected overlap set for [1000,2000]: {P_edgeStart, P_inside, P_edgeEnd} (3 peaks).
 */
public class RegulatoryPeakIndexTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        Path root = Files.createTempDirectory("regidx_test");
        try {
            writeFixture(root, "TESTE", "H3K27ac", new String[]{
                narrowPeakLine("chr1", 500,  900,  "P_before",    5.0),
                narrowPeakLine("chr1", 800,  1000, "P_edgeStart", 11.5),
                narrowPeakLine("chr1", 1200, 1300, "P_inside",    12.345),
                narrowPeakLine("chr1", 2000, 2100, "P_edgeEnd",   13.0),
                narrowPeakLine("chr1", 2200, 2300, "P_after",     14.0),
                narrowPeakLine("chr2", 1200, 1300, "P_otherChr",  15.0),
            });

            RegulatoryPeakIndex idx = new RegulatoryPeakIndex(root.toString());

            // ── Core overlap query, "chr1" form ──
            List<RegulatoryPeakIndex.Peak> hits = idx.peaksOverlapping("TESTE", "H3K27ac", "chr1", 1000, 2000);
            failures += check("window [1000,2000] returns exactly 3 peaks", hits.size() == 3);
            Set<Long> starts = new TreeSet<>();
            for (RegulatoryPeakIndex.Peak p : hits) starts.add(p.start);
            failures += check("window [1000,2000] returns the expected starts {800,1200,2000}",
                starts.equals(new TreeSet<>(Arrays.asList(800L, 1200L, 2000L))));
            failures += check("P_before (ends before window) excluded", !containsStart(hits, 500));
            failures += check("P_after (starts after window) excluded", !containsStart(hits, 2200));

            // ── Boundary cases explicitly ──
            failures += check("P_edgeStart (peak end == window start) is included",
                containsStart(hits, 800));
            failures += check("P_edgeEnd (peak start == window end) is included",
                containsStart(hits, 2000));

            // ── signalValue (real effect size) is preserved, not dropped ──
            RegulatoryPeakIndex.Peak inside = findStart(hits, 1200);
            failures += check("P_inside signalValue preserved", inside != null && inside.signalValue == 12.345);

            // ── chr/"chr" normalization: bare "1" must key the same bucket as "chr1" ──
            List<RegulatoryPeakIndex.Peak> hitsBare = idx.peaksOverlapping("TESTE", "H3K27ac", "1", 1000, 2000);
            failures += check("bare chromosome '1' normalizes the same as 'chr1'", hitsBare.size() == 3);

            // ── Different chromosome excluded — query a chr2 window that does NOT overlap
            // P_otherChr (chr2:1200-1300) to prove chr1's peaks never leak into a chr2 bucket. ──
            List<RegulatoryPeakIndex.Peak> hitsChr2 = idx.peaksOverlapping("TESTE", "H3K27ac", "chr2", 5000, 6000);
            failures += check("chr2 query outside its own peak returns empty (no chr1 leakage)", hitsChr2.isEmpty());
            // ── But a chr2 query overlapping P_otherChr's own span DOES return it, and only it. ──
            List<RegulatoryPeakIndex.Peak> hitsChr2Own = idx.peaksOverlapping("TESTE", "H3K27ac", "chr2", 1200, 1300);
            failures += check("chr2 query returns its own peak", hitsChr2Own.size() == 1);

            // ── Non-existent EID must return empty, not throw ──
            List<RegulatoryPeakIndex.Peak> noEid = idx.peaksOverlapping("NOPE", "H3K27ac", "chr1", 1000, 2000);
            failures += check("unknown EID returns empty (no exception)", noEid.isEmpty());

            // ── Non-existent mark for a known EID must return empty, not throw ──
            List<RegulatoryPeakIndex.Peak> noMark = idx.peaksOverlapping("TESTE", "H3K9ac", "chr1", 1000, 2000);
            failures += check("unknown mark returns empty (no exception)", noMark.isEmpty());

            // ── A query window entirely outside any peak returns empty ──
            List<RegulatoryPeakIndex.Peak> none = idx.peaksOverlapping("TESTE", "H3K27ac", "chr1", 5000, 6000);
            failures += check("non-overlapping window returns empty", none.isEmpty());

            // ── Corrupt file (not real gzip content) must not crash, just return empty ──
            Path badEidDir = root.resolve("BADEID");
            Files.createDirectories(badEidDir);
            Files.write(badEidDir.resolve("BADEID-H3K27ac.narrowPeak.gz"), "not actually gzip data".getBytes("UTF-8"));
            List<RegulatoryPeakIndex.Peak> corrupt = idx.peaksOverlapping("BADEID", "H3K27ac", "chr1", 1000, 2000);
            failures += check("corrupt/non-gzip file returns empty, does not throw", corrupt.isEmpty());

            if (failures == 0) {
                System.out.println("PASS: all RegulatoryPeakIndex tests passed");
            } else {
                System.out.println("FAIL: " + failures + " test(s) failed");
                System.exit(1);
            }
        } finally {
            deleteRecursive(root.toFile());
        }
    }

    private static String narrowPeakLine(String chr, long start, long end, String name, double signal) {
        // chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak — real narrowPeak
        // column order, matching the actual downloaded Roadmap Epigenomics files.
        return chr + "\t" + start + "\t" + end + "\t" + name + "\t500\t.\t" + signal + "\t20.0\t18.0\t50";
    }

    private static void writeFixture(Path root, String eid, String mark, String[] lines) throws IOException {
        Path dir = root.resolve(eid);
        Files.createDirectories(dir);
        File gz = dir.resolve(eid + "-" + mark + ".narrowPeak.gz").toFile();
        try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(gz))) {
            for (String line : lines) {
                out.write((line + "\n").getBytes("UTF-8"));
            }
        }
    }

    private static boolean containsStart(List<RegulatoryPeakIndex.Peak> peaks, long start) {
        return findStart(peaks, start) != null;
    }

    private static RegulatoryPeakIndex.Peak findStart(List<RegulatoryPeakIndex.Peak> peaks, long start) {
        for (RegulatoryPeakIndex.Peak p : peaks) if (p.start == start) return p;
        return null;
    }

    private static void deleteRecursive(File f) {
        File[] children = f.listFiles();
        if (children != null) for (File c : children) deleteRecursive(c);
        f.delete();
    }

    private static int check(String label, boolean cond) {
        if (!cond) {
            System.out.println("FAIL: " + label);
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
