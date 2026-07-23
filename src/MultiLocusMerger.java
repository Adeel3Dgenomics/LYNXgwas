import java.io.*;
import java.util.*;

/**
 * Locus Matrix phase 1: streams each selected dataset's raw GWAS file once,
 * keeping only the best (lowest) p-value seen per chr:pos position across
 * all datasets, filtered to p &lt;= pFilter. Writes the merged set to a
 * temp TSV (chrom, pos, p) for loci.LociIdentifier to consume unchanged.
 */
public class MultiLocusMerger {

    public static void merge(List<Config> configs, double pFilter, File outFile,
                              MultiLocusProgress progress) throws IOException {
        Map<String, double[]> best = new HashMap<>(); // "chr:pos" -> {p}

        progress.phase = "merge";
        progress.datasetTotal = configs.size();

        for (int i = 0; i < configs.size(); i++) {
            Config cfg = configs.get(i);
            progress.datasetIndex = i;
            progress.currentDataset = new File(cfg.outputDir).getName();
            streamOne(cfg, pFilter, best);
        }
        progress.datasetIndex = configs.size();

        outFile.getParentFile().mkdirs();
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile)))) {
            pw.println("chrom\tpos\tp");
            for (Map.Entry<String, double[]> e : best.entrySet()) {
                int c = e.getKey().indexOf(':');
                String chr = e.getKey().substring(0, c);
                String pos = e.getKey().substring(c + 1);
                pw.println(chr + "\t" + pos + "\t" + e.getValue()[0]);
            }
        }
    }

    private static void streamOne(Config cfg, double pFilter, Map<String, double[]> best) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(cfg.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) return;
            String[] cols = header.trim().split("\t");
            int iChr = GwasParser.colIdx(cols, cfg.colChr);
            int iPos = GwasParser.colIdx(cols, cfg.colPos);
            int iP   = GwasParser.colIdx(cols, cfg.colPvalue);
            if (iChr < 0 || iPos < 0 || iP < 0)
                throw new IOException("Missing required columns in " + cfg.gwasFile +
                    ": chr=" + cfg.colChr + " pos=" + cfg.colPos + " p=" + cfg.colPvalue);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = GwasParser.splitTab(line);
                if (f.length <= Math.max(iChr, Math.max(iPos, iP))) continue;

                double p;
                try { p = Double.parseDouble(f[iP].trim()); }
                catch (NumberFormatException e) { continue; }
                if (p <= 0 || p > pFilter) continue;

                String chr = f[iChr].trim();
                String pos = f[iPos].trim();
                String key = chr + ":" + pos;

                double[] cur = best.get(key);
                if (cur == null || p < cur[0]) {
                    best.put(key, new double[]{p});
                }
            }
        }
    }
}
