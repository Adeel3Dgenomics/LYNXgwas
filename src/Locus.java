import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

public class Locus {
    public final String id;
    public final int    index;
    public final String chr;
    public final long   start;
    public final long   end;
    public final long   paddedStart;
    public final long   paddedEnd;

    // Populated during GWAS streaming
    public final List<Snp> snps = new ArrayList<>();

    public Locus(int index, String chr, long start, long end, long padding) {
        this(UUID.randomUUID().toString(), index, chr, start, end, padding);
    }

    public Locus(String id, int index, String chr, long start, long end, long padding) {
        this.id          = id;
        this.index       = index;
        this.chr         = chr;
        this.start       = start;
        this.end         = end;
        this.paddedStart = Math.max(1, start - padding);
        this.paddedEnd   = end + padding;
    }

    public long mid() {
        return (start + end) / 2;
    }

    /** Numeric chromosome value for sorting (X=23, Y=24, MT=25). */
    public int chrInt() {
        return chrToInt(chr);
    }

    public static int chrToInt(String chr) {
        switch (chr.toUpperCase()) {
            case "X":  return 23;
            case "Y":  return 24;
            case "MT": case "M": return 25;
            default:
                try { return Integer.parseInt(chr); }
                catch (NumberFormatException e) { return 99; }
        }
    }

    @Override
    public String toString() {
        return "Locus " + index + " [" + id + " chr" + chr + ":" + start + "-" + end + "]";
    }
}
