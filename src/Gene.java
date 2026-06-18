import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Gene {
    public String geneName;
    public String geneId;
    public String chr;
    public String strand;
    public long   start;
    public long   end;
    public final List<Transcript> transcripts = new ArrayList<>();

    public Gene(String geneId, String geneName, String chr, String strand, long start, long end) {
        this.geneId   = geneId;
        this.geneName = geneName;
        this.chr      = chr;
        this.strand   = strand;
        this.start    = start;
        this.end      = end;
    }

    /**
     * Returns the "canonical" transcript: prefer protein_coding with the largest exon span,
     * else lncRNA with largest span.
     */
    public Transcript canonical() {
        Transcript best = transcripts.stream()
            .filter(t -> "protein_coding".equals(t.transcriptType))
            .max(Comparator.comparingLong(Transcript::exonSpan))
            .orElse(null);
        if (best == null) {
            best = transcripts.stream()
                .max(Comparator.comparingLong(Transcript::exonSpan))
                .orElse(null);
        }
        return best;
    }

    public long distanceTo(long pos) {
        if (pos < start) return start - pos;
        if (pos > end)   return pos - end;
        return 0;
    }
}
