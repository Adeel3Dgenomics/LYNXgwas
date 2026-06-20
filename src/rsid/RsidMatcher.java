package rsid;

import java.util.*;

// TODO: BUILD MISMATCH — verify/liftover build before production use.
// The local SNP database is hg19. Assumes GWAS data is also hg19.

public class RsidMatcher {

    /**
     * Pure function: match a single GWAS SNP against the position index.
     */
    public static MatchResult match(String chr, int pos, String ea, String nea,
                                    Map<Integer, List<DbSnpRecord>> posIndex) {
        String eaUp  = ea.toUpperCase();
        String neaUp = nea.toUpperCase();

        List<DbSnpRecord> candidates = posIndex.get(pos);
        if (candidates == null || candidates.isEmpty()) {
            return new MatchResult(chr, pos, ea, nea, null, "no_pos", null, 0);
        }

        List<DbSnpRecord> forwardMatches = new ArrayList<>();
        List<DbSnpRecord> reverseMatches = new ArrayList<>();

        for (DbSnpRecord rec : candidates) {
            // Forward: GWAS nea == ref AND GWAS ea in altList
            if (neaUp.equals(rec.ref) && rec.altList.contains(eaUp)) {
                forwardMatches.add(rec);
            }
            // Reverse: GWAS ea == ref AND GWAS nea in altList
            else if (eaUp.equals(rec.ref) && rec.altList.contains(neaUp)) {
                reverseMatches.add(rec);
            }
        }

        // TODO: optional strand-flip matching (A↔T, C↔G) — risks false matches on ambiguous SNPs

        if (!forwardMatches.isEmpty()) {
            DbSnpRecord best = selectNewestRsid(forwardMatches);
            return new MatchResult(chr, pos, ea, nea, best.rsid, "matched_forward",
                null, candidates.size());
        }
        if (!reverseMatches.isEmpty()) {
            DbSnpRecord best = selectNewestRsid(reverseMatches);
            return new MatchResult(chr, pos, ea, nea, best.rsid, "matched_reverse",
                null, candidates.size());
        }

        // Position matches but alleles don't — store best candidate rsID for diagnostics
        DbSnpRecord bestPosOnly = selectNewestRsid(candidates);
        return new MatchResult(chr, pos, ea, nea, null, "pos_only_allele_mismatch",
            bestPosOnly.rsid, candidates.size());
    }

    /**
     * Among matching candidates, pick the one from the newest dbSNP build.
     * If dbSnpBuild is absent from all candidates, pick the last one deterministically
     * (last record in VCF order at that position — typically the most recent addition).
     */
    static DbSnpRecord selectNewestRsid(List<DbSnpRecord> candidates) {
        if (candidates.size() == 1) return candidates.get(0);

        DbSnpRecord best = null;
        int highestBuild = -1;
        boolean anyHasBuild = false;

        for (DbSnpRecord rec : candidates) {
            if (rec.dbSnpBuild != null) {
                anyHasBuild = true;
                if (rec.dbSnpBuild > highestBuild) {
                    highestBuild = rec.dbSnpBuild;
                    best = rec;
                }
            }
        }

        // Fallback: if no build info, pick last record (most recent addition in VCF order)
        if (!anyHasBuild) {
            return candidates.get(candidates.size() - 1);
        }
        return best;
    }
}
