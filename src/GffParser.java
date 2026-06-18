import java.io.*;
import java.util.*;

/**
 * Streams the GFF3 (1.35 GB) once, building gene models for protein_coding
 * and lncRNA features only. Provides O(log n) interval-overlap queries.
 */
public class GffParser {

    /** Per-chromosome gene lists, sorted by start position for binary search. */
    private final Map<String, List<Gene>> chrGenes = new HashMap<>();
    private final Map<String, Transcript> txMap    = new HashMap<>();
    private final Map<String, Gene>       geneMap  = new HashMap<>();

    private static final Set<String> KEPT_TX_TYPES = new HashSet<>(
        Arrays.asList("protein_coding", "lncrna"));

    // -----------------------------------------------------------------------

    public static GffParser parse(Config config) throws IOException {
        GffParser gff = new GffParser();
        gff.stream(config.gff3File);
        return gff;
    }

    // -----------------------------------------------------------------------

    private void stream(String path) throws IOException {
        long t0 = System.currentTimeMillis();
        long lines = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(path), 4 * 1024 * 1024)) {
            String line;
            while ((line = br.readLine()) != null) {
                lines++;
                if (line.isEmpty() || line.charAt(0) == '#') continue;

                // Manually split on tab to avoid regex overhead
                int[] tabs = tabPositions(line);
                if (tabs.length < 8) continue; // need 8 tabs for 9 GFF3 columns

                // GFF3 cols: seqname(0) source(1) feature(2) start(3) end(4) score(5) strand(6) phase(7) attrs(8)
                // tabs[0..7] are the 8 tab characters between the 9 columns.
                String seqname = line.substring(0,           tabs[0]);  // col0
                String feature = line.substring(tabs[1] + 1, tabs[2]);  // col2
                long   fStart  = Long.parseLong(line.substring(tabs[2] + 1, tabs[3])); // col3
                long   fEnd    = Long.parseLong(line.substring(tabs[3] + 1, tabs[4])); // col4
                String strand  = line.substring(tabs[5] + 1, tabs[6]);  // col6
                String attrs   = line.substring(tabs[7] + 1);           // col8

                String chr = seqname.startsWith("chr") ? seqname.substring(3) : seqname;

                switch (feature) {
                    case "gene":
                        processGene(chr, strand, fStart, fEnd, attrs);
                        break;
                    case "transcript": case "mRNA":
                        processTranscript(strand, fStart, fEnd, attrs);
                        break;
                    case "exon":
                        addToTranscript(attrs, new Exon(fStart, fEnd), 0);
                        break;
                    case "CDS":
                        addToTranscript(attrs, new Exon(fStart, fEnd), 1);
                        break;
                    case "five_prime_UTR": case "three_prime_UTR": case "UTR":
                        addToTranscript(attrs, new Exon(fStart, fEnd, feature), 2);
                        break;
                }

                if (lines % 2_000_000 == 0) {
                    System.out.printf("[GffParser] %.0fM lines (%.1f s)%n",
                        lines / 1e6, (System.currentTimeMillis() - t0) / 1000.0);
                }
            }
        }

        // Sort each chromosome's gene list by start for overlap binary search
        for (List<Gene> genes : chrGenes.values()) {
            genes.sort(Comparator.comparingLong(g -> g.start));
        }
        txMap.clear(); // free memory — no longer needed

        int total = chrGenes.values().stream().mapToInt(List::size).sum();
        System.out.printf("[GffParser] Done: %.0fM lines, %d genes indexed (%.1f s)%n",
            lines / 1e6, total, (System.currentTimeMillis() - t0) / 1000.0);
    }

    // -----------------------------------------------------------------------

    private void processGene(String chr, String strand, long start, long end, String attrs) {
        Map<String, String> a = parseAttrs(attrs);
        String geneId   = a.getOrDefault("gene_id", a.getOrDefault("ID", ""));
        String geneName = a.getOrDefault("gene_name", geneId);
        if (geneId.isEmpty()) return;
        Gene gene = new Gene(geneId, geneName, chr, strand, start, end);
        geneMap.put(geneId, gene);
        chrGenes.computeIfAbsent(chr, k -> new ArrayList<>()).add(gene);
    }

    private void processTranscript(String strand, long start, long end, String attrs) {
        Map<String, String> a = parseAttrs(attrs);
        String txTypeRaw  = a.getOrDefault("transcript_type",
                            a.getOrDefault("transcript_biotype", ""));
        if (!KEPT_TX_TYPES.contains(txTypeRaw.toLowerCase())) return;

        String txId   = a.getOrDefault("transcript_id", a.getOrDefault("ID", ""));
        String parent = a.getOrDefault("Parent", "");
        if (txId.isEmpty()) return;

        Gene gene = findGene(parent);
        if (gene == null) return;

        Transcript tx = new Transcript(txId, txTypeRaw, strand, start, end);
        gene.transcripts.add(tx);
        txMap.put(txId, tx);
    }

    private void addToTranscript(String attrs, Exon exon, int kind) {
        Map<String, String> a = parseAttrs(attrs);
        String parent = a.getOrDefault("Parent", "");
        Transcript tx = txMap.get(parent);
        if (tx == null) return;
        switch (kind) {
            case 0: tx.exons.add(exon); break;
            case 1: tx.cds.add(exon);   break;
            case 2: tx.utrs.add(exon);  break;
        }
    }

    private Gene findGene(String id) {
        if (id.isEmpty()) return null;
        Gene g = geneMap.get(id);
        if (g != null) return g;
        int dot = id.lastIndexOf('.');
        return dot > 0 ? geneMap.get(id.substring(0, dot)) : null;
    }

    // -----------------------------------------------------------------------
    // Public query API

    /** All genes overlapping [queryStart, queryEnd] on the given chromosome. */
    public List<Gene> overlapping(String chr, long queryStart, long queryEnd) {
        List<Gene> genes = chrGenes.get(chr);
        if (genes == null || genes.isEmpty()) return Collections.emptyList();

        List<Gene> result = new ArrayList<>();
        int lo = 0, hi = genes.size() - 1;
        // Find the rightmost gene with start <= queryEnd
        int pivot = genes.size();
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            if (genes.get(mid).start <= queryEnd) { pivot = mid; lo = mid + 1; }
            else hi = mid - 1;
        }
        // Scan left from pivot to collect all overlapping
        for (int i = pivot; i >= 0; i--) {
            Gene g = genes.get(i);
            if (g.end < queryStart) break;
            if (g.end >= queryStart) result.add(g);
        }
        // Scan right from pivot+1 (genes with later starts that still overlap)
        for (int i = pivot + 1; i < genes.size(); i++) {
            Gene g = genes.get(i);
            if (g.start > queryEnd) break;
            result.add(g);
        }
        return result;
    }

    /** Names of up to 2 genes nearest to pos among candidates. */
    public List<String> nearestGeneNames(long pos, List<Gene> candidates) {
        if (candidates.isEmpty()) return Collections.emptyList();
        candidates.sort(Comparator.comparingLong(g -> g.distanceTo(pos)));
        List<String> names = new ArrayList<>();
        for (int i = 0; i < Math.min(2, candidates.size()); i++) {
            names.add(candidates.get(i).geneName);
        }
        return names;
    }

    // -----------------------------------------------------------------------
    // Helpers

    /** Returns indices of all tab characters in the line (up to 9). */
    private static int[] tabPositions(String line) {
        int[] pos = new int[9];
        int count = 0;
        for (int i = 0; i < line.length() && count < 9; i++) {
            if (line.charAt(i) == '\t') pos[count++] = i;
        }
        if (count < 9) return Arrays.copyOf(pos, count);
        return pos;
    }

    private static Map<String, String> parseAttrs(String attrs) {
        Map<String, String> map = new HashMap<>();
        int len = attrs.length(), i = 0;
        while (i < len) {
            int eq = attrs.indexOf('=', i);
            if (eq < 0) break;
            String key = attrs.substring(i, eq).trim();
            int semi = attrs.indexOf(';', eq + 1);
            String val = semi < 0 ? attrs.substring(eq + 1).trim()
                                  : attrs.substring(eq + 1, semi).trim();
            map.put(key, val);
            i = semi < 0 ? len : semi + 1;
        }
        return map;
    }
}
