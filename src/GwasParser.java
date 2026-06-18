import java.io.*;
import java.util.*;

/**
 * Streams the GWAS TSV once (sorted by chr:pos) and collects SNPs into each locus window.
 * Uses an active-loci approach: O(n + m) where n = GWAS rows, m = locus count.
 */
public class GwasParser {

    public static void parse(Config config, List<Locus> loci) throws IOException {
        if (loci.isEmpty()) return;

        // Sort loci by (chrInt, paddedStart) to match GWAS file order
        List<Locus> sorted = new ArrayList<>(loci);
        sorted.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.paddedStart));

        long t0 = System.currentTimeMillis();
        long lineCount = 0;
        long snpCount  = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(config.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) throw new IOException("GWAS file is empty");

            String[] cols = header.trim().split("\t");
            int iChr    = colIdx(cols, config.colChr);
            int iPos    = colIdx(cols, config.colPos);
            int iPval   = colIdx(cols, config.colPvalue);
            int iVarid  = config.colVarid.isEmpty() ? -1 : colIdx(cols, config.colVarid);
            int iRsid   = config.colRsid.isEmpty()  ? -1 : colIdx(cols, config.colRsid);
            int iEa     = colIdx(cols, config.colEa);
            int iNea    = colIdx(cols, config.colNea);

            if (iChr < 0 || iPos < 0 || iPval < 0) {
                throw new IOException("GWAS file missing required columns: "
                    + config.colChr + ", " + config.colPos + ", " + config.colPvalue);
            }

            // Active loci whose padded windows contain the current position
            List<Locus> activeLoci = new ArrayList<>();
            int nextLocusIdx = 0;         // pointer into sorted[]
            String currentChr = null;

            String line;
            while ((line = br.readLine()) != null) {
                lineCount++;
                if (line.isEmpty()) continue;

                String[] f = splitTab(line);
                if (f.length <= Math.max(iChr, Math.max(iPos, iPval))) continue;

                String chr = f[iChr].trim();
                long   pos;
                double pval;
                try {
                    pos  = Long.parseLong(f[iPos].trim());
                    pval = Double.parseDouble(f[iPval].trim());
                } catch (NumberFormatException e) {
                    continue; // skip malformed lines
                }

                // On chromosome change, flush active loci and reset
                if (!chr.equals(currentChr)) {
                    activeLoci.clear();
                    currentChr = chr;
                    // Advance nextLocusIdx to first locus on this chromosome
                    while (nextLocusIdx < sorted.size()
                           && sorted.get(nextLocusIdx).chrInt() < Locus.chrToInt(chr)) {
                        nextLocusIdx++;
                    }
                }

                // Add new loci whose padded windows begin at or before current pos
                while (nextLocusIdx < sorted.size()) {
                    Locus candidate = sorted.get(nextLocusIdx);
                    if (candidate.chrInt() != Locus.chrToInt(chr)) break;
                    if (candidate.paddedStart > pos) break;
                    activeLoci.add(candidate);
                    nextLocusIdx++;
                }

                // Remove loci whose padded windows ended before current pos
                activeLoci.removeIf(l -> l.paddedEnd < pos);

                if (activeLoci.isEmpty()) continue;

                // Build SNP object
                String varid = (iVarid >= 0 && iVarid < f.length) ? f[iVarid].trim() : chr + ":" + pos;
                String rsid  = (iRsid  >= 0 && iRsid  < f.length) ? f[iRsid].trim()  : "";
                String id    = (!rsid.isEmpty() && !rsid.equals(".")) ? rsid : varid;
                String ea    = (iEa  >= 0 && iEa  < f.length) ? f[iEa].trim()  : ".";
                String nea   = (iNea >= 0 && iNea < f.length) ? f[iNea].trim() : ".";

                Snp snp = new Snp(id, chr, pos, pval, ea, nea);

                for (Locus locus : activeLoci) {
                    locus.snps.add(snp);
                    snpCount++;
                }

                if (lineCount % 5_000_000 == 0) {
                    System.out.printf("[GwasParser] %.0fM lines processed, %d SNPs collected (%.1f s)%n",
                        lineCount / 1e6, snpCount, (System.currentTimeMillis() - t0) / 1000.0);
                }
            }
        }

        // Downsample: keep all p < 1e-4, random-sample remainder to maxSnpsPerLocus total
        int oversize = 0;
        Random rng = new Random(42);
        for (Locus locus : loci) {
            if (locus.snps.size() <= config.maxSnpsPerLocus) continue;
            List<Snp> significant = new ArrayList<>();
            List<Snp> rest        = new ArrayList<>();
            for (Snp s : locus.snps) {
                if (s.pvalue < 1e-4) significant.add(s);
                else                  rest.add(s);
            }
            int remaining = config.maxSnpsPerLocus - significant.size();
            if (remaining > 0 && !rest.isEmpty()) {
                Collections.shuffle(rest, rng);
                significant.addAll(rest.subList(0, Math.min(remaining, rest.size())));
            }
            locus.snps.clear();
            locus.snps.addAll(significant);
            oversize++;
        }

        System.out.printf("[GwasParser] Done: %.0fM lines, %d total SNPs collected (%.1f s)%s%n",
            lineCount / 1e6, snpCount, (System.currentTimeMillis() - t0) / 1000.0,
            oversize > 0 ? " [" + oversize + " loci downsampled]" : "");
    }

    // -----------------------------------------------------------------------

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++) {
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        }
        return -1;
    }

    /** Fast split on tab without regex overhead. */
    private static String[] splitTab(String line) {
        int count = 1;
        for (int i = 0; i < line.length(); i++) if (line.charAt(i) == '\t') count++;
        String[] parts = new String[count];
        int start = 0, idx = 0;
        for (int i = 0; i < line.length(); i++) {
            if (line.charAt(i) == '\t') {
                parts[idx++] = line.substring(start, i);
                start = i + 1;
            }
        }
        parts[idx] = line.substring(start);
        return parts;
    }
}
