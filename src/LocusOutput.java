import java.util.*;

/** Full output structure for one locus — serialised to locus_N.json. */
public class LocusOutput {
    public int    locusIndex;
    public String locusName;
    public String chr;
    public long   start;
    public long   end;
    public long   paddedStart;
    public long   paddedEnd;
    public String refPanel = "";    // population label, e.g. "EAS"

    public Snp          topSnp;
    public List<String> nearestGenes = new ArrayList<>();
    public List<Snp>    gwasSnps     = new ArrayList<>();
    public List<Gene>   genes        = new ArrayList<>();
    public LdTriangle   ldTriangle;   // null when LD unavailable

    public LocusContext locusContext;

    // ── LD triangle data ──────────────────────────────────────────────────────

    public static class LdTriangle {
        public int         snpCount;
        public int         cellSizePx  = 6;    // default triangle cell size in px
        public int         topSnpRank  = -1;   // rank (0-based index) of top SNP in snps list
        public List<LdSnp> snps   = new ArrayList<>();
        public double[][]  matrix;   // [snpCount][snpCount], upper triangle

        public static class LdSnp {
            public String  id;         // display ID (rsid if available)
            public String  bimId;      // ID from reference panel .bim file
            public long    pos;
            public boolean isTopSnp;
            public int     rank;       // position index in snps list (0-based)
        }
    }

    // ── Locus context ─────────────────────────────────────────────────────────

    public static class LocusContext {
        public LocusRef prevLocus;
        public LocusRef nextLocus;
    }

    public static class LocusRef {
        public int    index;
        public String chr;
        public long   start;
        public long   end;
        public long   mid;
        public long   distanceBp;

        public LocusRef(int index, String chr, long start, long end, long mid, long distanceBp) {
            this.index      = index;
            this.chr        = chr;
            this.start      = start;
            this.end        = end;
            this.mid        = mid;
            this.distanceBp = distanceBp;
        }
    }
}
