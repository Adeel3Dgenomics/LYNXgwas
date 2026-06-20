public class Snp {
    public String id;          // rsid if available, else varid (chr:pos:ref:alt)
    public String bimId;       // ID as it appears in the reference panel .bim file (may differ from id)
    public String chr;
    public long   pos;
    public double pvalue;
    public double negLog10P;
    public double r2 = -1;     // r² with top SNP; -1 = no LD data
    public String ea;
    public String nea;

    // Optional GWAS columns — NaN means not configured/available
    public double beta      = Double.NaN;
    public double oddsRatio = Double.NaN;
    public double se        = Double.NaN;
    public double sampleN   = Double.NaN;
    public double maf       = Double.NaN;
    public double infoScore = Double.NaN;

    public Snp(String id, String chr, long pos, double pvalue, String ea, String nea) {
        this.id        = id;
        this.chr       = chr;
        this.pos       = pos;
        this.pvalue    = pvalue;
        this.negLog10P = pvalue > 0 ? -Math.log10(pvalue) : 308.0;
        this.ea        = ea;
        this.nea       = nea;
    }
}
