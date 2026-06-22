
import java.util.*;

/**
 * Stable SNP identifier: chr:pos:ref:alt (alleles alphabetically sorted).
 * Every analysis artifact keys SNPs by this format so results
 * are joinable across tools regardless of allele order conventions.
 */
public class StableSnpId {

    public final String chr;
    public final long pos;
    public final String a1; // alphabetically first
    public final String a2; // alphabetically second
    public final String id;

    private StableSnpId(String chr, long pos, String a1, String a2) {
        this.chr = chr;
        this.pos = pos;
        this.a1 = a1;
        this.a2 = a2;
        this.id = chr + ":" + pos + ":" + a1 + ":" + a2;
    }

    public static StableSnpId of(String chr, long pos, String alleleA, String alleleB) {
        String a = alleleA.toUpperCase();
        String b = alleleB.toUpperCase();
        return a.compareTo(b) <= 0
            ? new StableSnpId(chr, pos, a, b)
            : new StableSnpId(chr, pos, b, a);
    }

    public static StableSnpId parse(String id) {
        String[] parts = id.split(":");
        if (parts.length < 4) return null;
        try {
            return new StableSnpId(parts[0], Long.parseLong(parts[1]), parts[2], parts[3]);
        } catch (NumberFormatException e) { return null; }
    }

    @Override public String toString() { return id; }
    @Override public int hashCode() { return id.hashCode(); }
    @Override public boolean equals(Object o) {
        return o instanceof StableSnpId && id.equals(((StableSnpId) o).id);
    }
}
