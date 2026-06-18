import java.util.ArrayList;
import java.util.List;

public class Transcript {
    public String transcriptId;
    public String transcriptType;
    public String strand;
    public long   start;
    public long   end;
    public final List<Exon> exons = new ArrayList<>();
    public final List<Exon> cds   = new ArrayList<>();
    public final List<Exon> utrs  = new ArrayList<>();

    public Transcript(String transcriptId, String transcriptType, String strand, long start, long end) {
        this.transcriptId   = transcriptId;
        this.transcriptType = transcriptType;
        this.strand         = strand;
        this.start          = start;
        this.end            = end;
    }

    /** Total exon span — used to pick the "best" transcript per gene. */
    public long exonSpan() {
        return exons.stream().mapToLong(e -> e.end - e.start).sum();
    }
}
