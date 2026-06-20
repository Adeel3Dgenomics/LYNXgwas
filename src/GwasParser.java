import java.io.*;
import java.util.*;

/**
 * Streams the GWAS TSV once (sorted by chr:pos) and collects SNPs into each locus window.
 * Uses an active-loci approach: O(n + m) where n = GWAS rows, m = locus count.
 */
public class GwasParser {

    public static void parse(Config config, List<Locus> loci) throws IOException {
        if (loci.isEmpty()) return;

        List<Locus> sorted = new ArrayList<>(loci);
        sorted.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.paddedStart));

        long t0 = System.currentTimeMillis();
        long lineCount = 0;
        long snpCount  = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(config.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) throw new IOException("GWAS file is empty");

            String[] cols = header.trim().split("\t");
            String avail = String.join(", ", cols);

            // Required columns — fail fast with clear error if not found
            int iChr  = requireCol(cols, "col.chr",    config.colChr,    "chromosome", avail);
            int iPos  = requireCol(cols, "col.pos",    config.colPos,    "position",   avail);
            int iPval = requireCol(cols, "col.pvalue", config.colPvalue, "p-value",    avail);
            int iEa   = requireCol(cols, "col.ea",     config.colEa,    "effect allele (A1)", avail);
            int iNea  = requireCol(cols, "col.nea",    config.colNea,   "other allele (A2)",  avail);

            // Optional identity columns
            int iVarid = optionalCol(cols, "col.varid", config.colVarid, "variant ID", avail);
            int iRsid  = optionalCol(cols, "col.rsid",  config.colRsid,  "rsID",       avail);

            // Optional GWAS statistic columns
            int iBeta = optionalCol(cols, "col.beta", config.colBeta, "beta",             avail);
            int iOr   = optionalCol(cols, "col.or",   config.colOr,   "odds ratio",       avail);
            int iSe   = optionalCol(cols, "col.se",   config.colSe,   "standard error",   avail);
            int iN    = optionalCol(cols, "col.n",    config.colN,    "sample size",       avail);
            int iMaf  = optionalCol(cols, "col.maf",  config.colMaf,  "minor allele freq", avail);
            int iInfo = optionalCol(cols, "col.info", config.colInfo, "imputation info",   avail);

            List<Locus> activeLoci = new ArrayList<>();
            int nextLocusIdx = 0;
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
                    continue;
                }

                if (!chr.equals(currentChr)) {
                    activeLoci.clear();
                    currentChr = chr;
                    while (nextLocusIdx < sorted.size()
                           && sorted.get(nextLocusIdx).chrInt() < Locus.chrToInt(chr)) {
                        nextLocusIdx++;
                    }
                }

                while (nextLocusIdx < sorted.size()) {
                    Locus candidate = sorted.get(nextLocusIdx);
                    if (candidate.chrInt() != Locus.chrToInt(chr)) break;
                    if (candidate.paddedStart > pos) break;
                    activeLoci.add(candidate);
                    nextLocusIdx++;
                }

                activeLoci.removeIf(l -> l.paddedEnd < pos);
                if (activeLoci.isEmpty()) continue;

                // Build SNP
                String varid = (iVarid >= 0 && iVarid < f.length) ? f[iVarid].trim() : chr + ":" + pos;
                String rsid  = (iRsid  >= 0 && iRsid  < f.length) ? f[iRsid].trim()  : "";
                String id    = (!rsid.isEmpty() && !rsid.equals(".")) ? rsid : varid;
                String ea    = (iEa  < f.length) ? f[iEa].trim()  : ".";
                String nea   = (iNea < f.length) ? f[iNea].trim() : ".";

                Snp snp = new Snp(id, chr, pos, pval, ea, nea);

                // Optional statistic columns
                if (iBeta >= 0) snp.beta      = parseOptionalDouble(f, iBeta);
                if (iOr   >= 0) snp.oddsRatio = parseOptionalDouble(f, iOr);
                if (iSe   >= 0) snp.se        = parseOptionalDouble(f, iSe);
                if (iN    >= 0) snp.sampleN   = parseOptionalDouble(f, iN);
                if (iMaf  >= 0) snp.maf       = parseOptionalDouble(f, iMaf);
                if (iInfo >= 0) snp.infoScore  = parseOptionalDouble(f, iInfo);

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

    // ── Column validation ─────────────────────────────────────────────────

    private static int requireCol(String[] header, String configKey, String colName,
                                  String description, String available) throws IOException {
        int idx = colIdx(header, colName);
        if (idx < 0) {
            throw new IOException(String.format(
                "Required GWAS column not found: %s (configured as %s='%s'). " +
                "Available columns: [%s]", description, configKey, colName, available));
        }
        return idx;
    }

    private static int optionalCol(String[] header, String configKey, String colName,
                                   String description, String available) throws IOException {
        if (colName == null || colName.isEmpty()) return -1;
        int idx = colIdx(header, colName);
        if (idx < 0) {
            throw new IOException(String.format(
                "Optional GWAS column configured but not found in header: %s " +
                "(configured as %s='%s'). Available columns: [%s]. " +
                "Remove the mapping or correct the column name.",
                description, configKey, colName, available));
        }
        return idx;
    }

    static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++) {
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        }
        return -1;
    }

    private static double parseOptionalDouble(String[] fields, int idx) {
        if (idx >= fields.length) return Double.NaN;
        String v = fields[idx].trim();
        if (v.isEmpty() || v.equals(".") || v.equalsIgnoreCase("NA") || v.equalsIgnoreCase("nan"))
            return Double.NaN;
        try { return Double.parseDouble(v); }
        catch (NumberFormatException e) { return Double.NaN; }
    }

    /** Fast split on tab without regex overhead. */
    static String[] splitTab(String line) {
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
