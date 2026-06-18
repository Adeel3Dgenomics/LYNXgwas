import java.io.*;
import java.util.*;

/**
 * Precomputes a fixed-size binned genome-wide Manhattan skyline.
 *
 * Streams the GWAS file once, assigns each SNP to a genomic bin using
 * GRCh38 chromosome lengths, and stores max(-log10 P) per bin.
 * Output: a small JSON array (~3000 entries) written to genome_skyline.json/js.
 *
 * The chromosome offset table is written into the JSON so the front-end
 * uses the exact same coordinate system for locus highlights.
 */
public class GenomeSkyline {

    static final int DEFAULT_BIN_COUNT = 3000;

    // GRCh38 chromosome lengths (bp)
    static final String[] CHR_ORDER = {
        "1","2","3","4","5","6","7","8","9","10",
        "11","12","13","14","15","16","17","18","19","20",
        "21","22","X"
    };
    static final long[] CHR_LENGTHS = {
        248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
        170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
        135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
         90338345L,  83257441L,  80373285L,  58617616L,  64444167L,
         46709983L,  50818468L, 156040895L
    };

    // Cumulative offsets: chrOffset[i] = sum of lengths of chr 0..i-1
    static final long[] CHR_OFFSETS = new long[CHR_ORDER.length];
    static final long GENOME_LENGTH;
    static {
        long cum = 0;
        for (int i = 0; i < CHR_ORDER.length; i++) {
            CHR_OFFSETS[i] = cum;
            cum += CHR_LENGTHS[i];
        }
        GENOME_LENGTH = cum;
    }

    static int chrIndex(String chr) {
        String c = chr.replaceFirst("(?i)^chr", "");
        for (int i = 0; i < CHR_ORDER.length; i++)
            if (CHR_ORDER[i].equalsIgnoreCase(c)) return i;
        return -1;
    }

    public static void generate(Config config) throws IOException {
        int binCount = DEFAULT_BIN_COUNT;
        double binSize = (double) GENOME_LENGTH / binCount;

        double[] bins = new double[binCount]; // max -log10(P) per bin
        long snpCount = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(config.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) return;

            String[] cols = header.trim().split("\t");
            int iChr  = colIdx(cols, config.colChr);
            int iPos  = colIdx(cols, config.colPos);
            int iPval = colIdx(cols, config.colPvalue);
            if (iChr < 0 || iPos < 0 || iPval < 0) {
                System.err.println("[GenomeSkyline] Missing required columns — skipping");
                return;
            }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = splitTab(line);
                if (f.length <= Math.max(iChr, Math.max(iPos, iPval))) continue;

                String chr = f[iChr].trim();
                int ci = chrIndex(chr);
                if (ci < 0) continue;

                long pos;
                double pval;
                try {
                    pos  = Long.parseLong(f[iPos].trim());
                    pval = Double.parseDouble(f[iPval].trim());
                } catch (NumberFormatException e) { continue; }

                if (pval <= 0 || Double.isNaN(pval)) continue;

                double nlp = -Math.log10(pval);
                long cumPos = CHR_OFFSETS[ci] + pos;
                int bi = (int) Math.min(binCount - 1, cumPos / (long) Math.ceil(binSize));
                if (nlp > bins[bi]) bins[bi] = nlp;
                snpCount++;
            }
        }

        // Write JSON
        String dataDir = config.outputDir + "/data";
        new File(dataDir).mkdirs();

        StringBuilder sb = new StringBuilder(binCount * 20 + 512);
        sb.append("{\"bin_count\":").append(binCount);
        sb.append(",\"genome_length\":").append(GENOME_LENGTH);
        sb.append(",\"snp_count\":").append(snpCount);

        // Chromosome offset table
        sb.append(",\"chromosomes\":[");
        for (int i = 0; i < CHR_ORDER.length; i++) {
            if (i > 0) sb.append(',');
            sb.append("{\"chr\":\"").append(CHR_ORDER[i]).append("\"");
            sb.append(",\"offset\":").append(CHR_OFFSETS[i]);
            sb.append(",\"length\":").append(CHR_LENGTHS[i]).append('}');
        }
        sb.append(']');

        // Bins: only emit non-zero bins to save space
        sb.append(",\"bins\":[");
        boolean first = true;
        for (int i = 0; i < binCount; i++) {
            if (bins[i] <= 0) continue;
            if (!first) sb.append(',');
            first = false;
            long binStart = (long) (i * binSize);
            sb.append('[').append(i).append(',');
            sb.append(String.format("%.2f", bins[i])).append(']');
        }
        sb.append("]}");

        String json = sb.toString();

        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/genome_skyline.json")))) {
            pw.print(json);
        }
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(dataDir + "/genome_skyline.js")))) {
            pw.print("(function(){window.GENOME_SKYLINE=" + json + ";})();");
        }

        System.out.printf("[GenomeSkyline] %d bins, %d SNPs, payload %.1f KB%n",
            binCount, snpCount, json.length() / 1024.0);
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }

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
