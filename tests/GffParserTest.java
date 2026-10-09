import java.io.*;
import java.util.*;

/**
 * Standalone regression test (no external test framework) for GffParser.overlapping.
 *
 * The old search scanned left from the last gene starting before the query and stopped at the first gene ending
 * before it. Genes are sorted by start, not end, so a long gene starting further left that still covers the query
 * was missed whenever a short gene sat in between (29% of loci lost genes: MAGMA, gene tracks, nearest genes).
 */
public class GffParserTest {

    public static void main(String[] args) throws Exception {
        File f = File.createTempFile("gff_test", ".gff3");
        try (PrintWriter pw = new PrintWriter(f)) {
            pw.println("##gff-version 3");
            pw.println("chr1\tT\tgene\t1000\t100000\t.\t+\t.\tID=G1;gene_id=G1;gene_name=LONG");
            pw.println("chr1\tT\tgene\t2000\t3000\t.\t+\t.\tID=G2;gene_id=G2;gene_name=SHORT");
            pw.println("chr1\tT\tgene\t200000\t210000\t.\t-\t.\tID=G3;gene_id=G3;gene_name=NEXT");
        }
        Config c = new Config();
        c.gff3File = f.getAbsolutePath();
        GffParser gff = GffParser.parse(c);
        f.delete();

        int failures = 0;
        failures += check("query inside a long gene, short gene in between: the long gene is found",
            names(gff.overlapping("1", 50000, 60000)), "LONG");
        failures += check("query covering both genes: both", names(gff.overlapping("1", 2500, 2600)), "LONG SHORT");
        failures += check("query between genes: none", names(gff.overlapping("1", 150000, 160000)), "");
        failures += check("query before every gene: none, and no exception", names(gff.overlapping("1", 10, 20)), "");
        failures += check("query after every gene: none", names(gff.overlapping("1", 300000, 400000)), "");
        failures += check("query spanning all: all three", names(gff.overlapping("1", 0, 1_000_000)), "LONG NEXT SHORT");
        failures += check("unknown chromosome: none", names(gff.overlapping("2", 0, 1_000_000)), "");

        if (failures == 0) System.out.println("PASS: all GffParser tests passed");
        else { System.out.println("FAIL: " + failures + " test(s) failed"); System.exit(1); }
    }

    /** Gene names, sorted (the test is about which genes, not their order). */
    private static String names(List<Gene> genes) {
        List<String> n = new ArrayList<>();
        for (Gene g : genes) n.add(g.geneName);
        Collections.sort(n);
        return String.join(" ", n);
    }

    private static int check(String label, String actual, String expected) {
        if (!actual.equals(expected)) {
            System.out.println("FAIL: " + label + " — expected [" + expected + "], got [" + actual + "]");
            return 1;
        }
        System.out.println("PASS: " + label);
        return 0;
    }
}
