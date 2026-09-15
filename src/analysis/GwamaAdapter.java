
import java.io.*;
import java.util.*;

/**
 * GWAMA meta-analysis adapter.
 *
 * GWAMA (Magi &amp; Morris 2010) meta-analyzes two or more per-cohort summary
 * statistics files. LYNXgwas currently drives this per single project/locus, so this
 * adapter writes exactly one cohort file in GWAMA's standard quantitative-trait
 * column layout (MARKERNAME/EA/NEA/BETA/SE/N/EAF/STRAND) from the locus's own
 * harmonized GWAS. Run on a single cohort this is a pass-through (genomic-control
 * recalibration of one dataset, not a true meta-analysis) — it exists so the tool
 * is functional end-to-end rather than failing on a missing input file; real
 * multi-cohort combination needs a cross-project runner (see SusiexAdapter for the
 * equivalent cross-project pattern used for SuSiEx).
 */
public class GwamaAdapter {

    public static void prepareRun(File harmonizedDir, File runDir) throws IOException {
        runDir.mkdirs();
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists()) throw new IOException("harmonized_gwas.tsv not found. Run the base pipeline first.");

        File cohortFile = new File(runDir, "cohort1.txt");
        File listFile = new File(runDir, "gwama_input.txt");

        int count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(cohortFile)))) {
            pw.println("MARKERNAME\tEA\tNEA\tBETA\tSE\tN\tEAF\tSTRAND");
            br.readLine(); // header: snp_id chr pos ea nea pvalue beta se or n maf info rsid varid
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;
                String snp = f[0], ea = f[3], nea = f[4], beta = f[6], se = f[7];
                if (beta.equals("NA") || se.equals("NA")) continue;
                String n = (f.length > 9 && !f[9].equals("NA") && !f[9].isEmpty()) ? f[9] : "0";
                String eaf = (f.length > 10 && !f[10].equals("NA")) ? f[10] : "NA";
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t+%n", snp, ea, nea, beta, se, n, eaf);
                count++;
            }
        }
        if (count == 0) throw new IOException("No usable SNPs (beta/se) in harmonized_gwas.tsv for GWAMA input.");

        // GWAMA's -i flag takes a list of cohort file paths, one per line.
        try (PrintWriter pw = new PrintWriter(new FileWriter(listFile))) {
            pw.println(cohortFile.getAbsolutePath().replace("\\", "/"));
        }

        System.out.printf("[GwamaAdapter] Wrote %d SNPs to single-cohort GWAMA input (%s)%n",
            count, cohortFile.getName());
    }
}
