package rsid;

public class MatchResult {
    public final String chr;
    public final int pos;
    public final String gwasEa;
    public final String gwasNea;
    public final String assignedRsid;
    public final String matchReason;
    public final String posOnlyRsid;
    public final int nCandidatesAtPos;

    public MatchResult(String chr, int pos, String gwasEa, String gwasNea,
                       String assignedRsid, String matchReason,
                       String posOnlyRsid, int nCandidatesAtPos) {
        this.chr              = chr;
        this.pos              = pos;
        this.gwasEa           = gwasEa;
        this.gwasNea          = gwasNea;
        this.assignedRsid     = assignedRsid;
        this.matchReason      = matchReason;
        this.posOnlyRsid      = posOnlyRsid;
        this.nCandidatesAtPos = nCandidatesAtPos;
    }

    public String toCsv() {
        return String.join(",", chr, String.valueOf(pos), gwasEa, gwasNea,
            assignedRsid != null ? assignedRsid : "",
            matchReason,
            posOnlyRsid != null ? posOnlyRsid : "",
            String.valueOf(nCandidatesAtPos));
    }

    public static String csvHeader() {
        return "chr,pos,gwas_ea,gwas_nea,assigned_rsid,match_reason,pos_only_rsid,n_candidates_at_pos";
    }
}
