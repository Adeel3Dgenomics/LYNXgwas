public class Exon {
    public final long start;
    public final long end;
    public final String type; // null for exon/CDS; "five_prime_UTR" or "three_prime_UTR" for UTRs

    public Exon(long start, long end) {
        this.start = start;
        this.end   = end;
        this.type  = null;
    }

    public Exon(long start, long end, String type) {
        this.start = start;
        this.end   = end;
        this.type  = type;
    }
}
