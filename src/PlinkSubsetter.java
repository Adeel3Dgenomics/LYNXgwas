import java.io.*;
import java.util.*;

/**
 * Phase 2: sequentially extract one PLINK bfile subset per locus.
 * Builds a chr:pos → bimId map from each .bim file so LdCalculator
 * can match GWAS SNPs to reference-panel variants without relying on rsid.
 */
public class PlinkSubsetter {

    public static class SubsetResult {
        /** "chr:pos" → BIM variant ID (column 2 of .bim). */
        public final Map<String, String> chrPosToVarid = new LinkedHashMap<>();
        /** BIM ID of the top SNP; null if top SNP not found in reference panel. */
        public String topSnpBimId = null;
        /** True if PLINK extraction succeeded and .bim exists. */
        public boolean ok = false;
    }

    public static Map<Integer, SubsetResult> subsetAll(
            List<Locus>         loci,
            Map<Integer, Snp>   topSnps,
            String              plinkBin,
            Config              config,
            ProgressTracker     progress) throws Exception {

        Map<Integer, SubsetResult> results = new LinkedHashMap<>();
        new File(config.plinkSubsetsDir()).mkdirs();

        for (Locus locus : loci) {
            if (progress != null)
                progress.update("Extracting PLINK subsets", locus.index, loci.size());

            SubsetResult result = subsetLocus(locus, plinkBin, config);

            // Identify top SNP BIM ID
            Snp topSnp = topSnps.get(locus.index);
            if (topSnp != null && result.ok) {
                String key = topSnp.chr + ":" + topSnp.pos;
                result.topSnpBimId = result.chrPosToVarid.get(key);
                if (result.topSnpBimId == null) {
                    System.err.printf("[WARN] Locus %d: top SNP %s (%s) absent from ref panel%n",
                        locus.index, topSnp.id, key);
                } else {
                    topSnp.bimId = result.topSnpBimId;
                }
            }

            results.put(locus.index, result);
            if (result.ok)
                System.out.printf("[Subset] Locus %d: %d variants extracted%n",
                    locus.index, result.chrPosToVarid.size());
        }
        return results;
    }

    private static SubsetResult subsetLocus(Locus locus, String plinkBin, Config config)
            throws IOException, InterruptedException {

        SubsetResult result = new SubsetResult();
        String prefix = config.plinkSubsetsDir() + "/locus_" + locus.index;

        List<String> cmd = Arrays.asList(
            plinkBin,
            "--bfile",   config.refPanelPath,
            "--chr",     locus.chr,
            "--from-bp", String.valueOf(locus.paddedStart),
            "--to-bp",   String.valueOf(locus.paddedEnd),
            "--make-bed",
            "--out",     prefix,
            "--silent"
        );

        Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
            while (br.readLine() != null) { /* drain */ }
        }
        int exit = proc.waitFor();

        File bimFile = new File(prefix + ".bim");
        if (exit != 0 || !bimFile.exists()) {
            System.err.printf("[WARN] Locus %d: PLINK subset failed (exit=%d)%n", locus.index, exit);
            return result;
        }

        // Parse .bim: chr  snpId  morgans  pos  a1  a2  (tab-delimited)
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 4) continue;
                String chr   = f[0].replaceFirst("^chr", "");
                String bimId = f[1];
                String pos   = f[3].trim();
                if (!pos.isEmpty() && !pos.equals("0"))
                    result.chrPosToVarid.put(chr + ":" + pos, bimId);
            }
        }

        // Delete .log and .nosex to save space; keep .bed/.bim/.fam for LD
        new File(prefix + ".log").delete();
        new File(prefix + ".nosex").delete();
        result.ok = true;
        return result;
    }

    /** Find PLINK 1.9 executable. */
    public static String findPlink(Config config) {
        String[] candidates = {"bin/plink.exe", "bin/plink", "plink.exe", "plink"};
        for (String c : candidates) {
            if (new File(c).exists()) return new File(c).getAbsolutePath();
        }
        return null;
    }
}
