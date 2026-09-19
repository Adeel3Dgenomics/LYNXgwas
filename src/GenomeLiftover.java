import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import java.io.File;

/**
 * hg19 (GRCh37) -> hg38 (GRCh38) coordinate conversion, needed because the GWAS Catalog's bulk
 * associations file is GRCh38-only (verified: EBI publishes no GRCh37 variant of it, and rs2238057's
 * position matches Open Targets' explicitly-GRCh38 coordinate exactly) while LYNXgwas projects can
 * legitimately be configured as either build — GRCh37 is in fact the project's own default. Without
 * this, a GRCh37 project's known/novel check would silently compare hg19 locus coordinates against
 * hg38 Catalog positions with no warning at all.
 *
 * Uses HTSJDK's LiftOver (already a dependency for VCF/dbSNP work — no new dependency here) against the
 * standard UCSC hg19ToHg38.over.chain file. Validated against two independently-known anchors before
 * wiring this in: rs429358's well-documented hg19 position (chr19:45,411,941) lifts to exactly
 * chr19:44,908,684 — the same hg38 position Open Targets independently reports — and rs2238057's known
 * hg38 position round-trips through hg38->hg19->hg38 back to the exact original coordinate.
 *
 * L2G/enhancer-to-gene lookups do NOT need this: they resolve by rsID through Open Targets' own
 * search(), which is build-agnostic. Only the GWAS-Catalog region-overlap check compares raw chr:pos
 * numbers directly, so that's the only place a build mismatch can silently corrupt results.
 */
public class GenomeLiftover {

    private static final String CHAIN_FILE = "resources/hg19ToHg38.over.chain";
    private static volatile LiftOver liftOver;
    private static volatile String initError;

    private static LiftOver get() {
        if (liftOver != null || initError != null) return liftOver;
        synchronized (GenomeLiftover.class) {
            if (liftOver != null || initError != null) return liftOver;
            File f = new File(CHAIN_FILE);
            if (!f.isFile()) {
                initError = "hg19ToHg38 chain file not found at " + CHAIN_FILE +
                    " — download it from hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz and unzip it there.";
                System.err.println("[GenomeLiftover] " + initError);
                return null;
            }
            try {
                liftOver = new LiftOver(f);
            } catch (Exception e) {
                initError = "Failed to load liftover chain: " + e.getMessage();
                System.err.println("[GenomeLiftover] " + initError);
            }
        }
        return liftOver;
    }

    public static class Result {
        public boolean ok;
        public String chr;
        public long start, end;
        public String error;
    }

    /** True if this build string means the coordinates need converting before comparing against the
     *  (GRCh38-only) GWAS Catalog index. */
    public static boolean needsLiftover(String genomeBuild) {
        return genomeBuild != null && genomeBuild.trim().equalsIgnoreCase("GRCh37");
    }

    /** Converts a locus region from GRCh37 to GRCh38. chr may be with or without a "chr" prefix.
     *  result.ok is false (with .error set) if the chain isn't loaded or either endpoint has no
     *  confident mapping in this region — callers should treat that as "cannot verify", not silently
     *  fall back to the un-lifted coordinates. */
    public static Result toGRCh38(String chr, long start, long end) {
        Result r = new Result();
        LiftOver lo = get();
        if (lo == null) {
            r.error = initError;
            return r;
        }
        String ucscChr = chr.regionMatches(true, 0, "chr", 0, 3) ? chr : "chr" + chr;
        Interval lifted;
        try {
            lifted = lo.liftOver(new Interval(ucscChr, (int) start, (int) end));
        } catch (Exception e) {
            r.error = "Liftover threw: " + e.getMessage();
            return r;
        }
        if (lifted == null) {
            r.error = "No confident GRCh37->GRCh38 mapping for " + chr + ":" + start + "-" + end;
            return r;
        }
        r.ok = true;
        r.chr = lifted.getContig().replaceFirst("(?i)^chr", "");
        r.start = lifted.getStart();
        r.end = lifted.getEnd();
        return r;
    }
}
