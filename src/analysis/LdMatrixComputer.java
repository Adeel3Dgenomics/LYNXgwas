
import java.io.*;
import java.util.*;

/**
 * Step 4: Compute LD matrices (r, r², D') from the harmonized/matched panel,
 * in harmonized effect-allele orientation, with an authoritative ld_snp_order.txt.
 *
 * Uses PLINK 1.9 for r and r² computation.
 * D' requires a separate PLINK --ld-window command or post-processing.
 */
public class LdMatrixComputer {

    public static class Result {
        public boolean ok;
        public String error;
        public int snpCount;
        public String inputHash;
    }

    public static String computeInputHash(File matchedDir) throws IOException {
        return ContentHasher.combineHashes(
            ContentHasher.hashFile(new File(matchedDir, "matched_ref.bim")),
            ContentHasher.hashFile(new File(matchedDir, "matched_ref.bed")),
            ContentHasher.hashFile(new File(matchedDir, "matched_ref.fam"))
        );
    }

    public static final int DEFAULT_LD_WINDOW = 200;

    public static Result run(File matchedDir, File harmonizedDir, File ldDir,
                              Config config) throws IOException {
        return run(matchedDir, harmonizedDir, ldDir, config, DEFAULT_LD_WINDOW);
    }

    public static Result run(File matchedDir, File harmonizedDir, File ldDir,
                              Config config, int ldWindow) throws IOException {
        Result result = new Result();
        ldDir.mkdirs();

        String inputHash = ContentHasher.combineHashes(
            computeInputHash(matchedDir), String.valueOf(ldWindow));
        result.inputHash = inputHash;

        StepManifest existing = StepManifest.read(ldDir, "ld");
        if (existing != null && existing.isCurrent(inputHash)) {
            if (new File(ldDir, "ld_r2.matrix").exists() &&
                new File(ldDir, "ld_snp_order.txt").exists()) {
                result.ok = true;
                return result;
            }
        }

        String plinkBin = PlinkSubsetter.findPlink(config);
        if (plinkBin == null) {
            result.error = "PLINK binary not found";
            return result;
        }

        File refBed = new File(matchedDir, "matched_ref.bed");
        if (!refBed.exists()) {
            result.error = "matched_ref.bed not found";
            return result;
        }

        // ── Read all BIM entries ──
        File bimFile = new File(matchedDir, "matched_ref.bim");
        List<String> allBimIds = new ArrayList<>();
        List<String> allSnpKeys = new ArrayList<>();
        List<Long> allPositions = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;
                String chr = f[0].replaceFirst("^chr", "");
                String bimId = f[1];
                long pos = Long.parseLong(f[3].trim());
                String a1 = f[4].trim().toUpperCase();
                String a2 = f[5].trim().toUpperCase();
                allBimIds.add(bimId);
                allSnpKeys.add(chr + ":" + pos + ":" + a1 + ":" + a2);
                allPositions.add(pos);
            }
        }

        if (allBimIds.isEmpty()) {
            result.error = "No SNPs in matched ref panel";
            return result;
        }

        // ── Find lead SNP (lowest p-value) from harmonized GWAS ──
        int leadIdx = findLeadSnpIndex(harmonizedDir, allPositions);
        if (leadIdx < 0) leadIdx = allBimIds.size() / 2;

        // ── Select window around lead SNP ──
        int winStart = Math.max(0, leadIdx - ldWindow);
        int winEnd = Math.min(allBimIds.size(), leadIdx + ldWindow + 1);
        List<String> windowBimIds = allBimIds.subList(winStart, winEnd);
        List<String> windowSnpKeys = allSnpKeys.subList(winStart, winEnd);

        System.out.printf("[LdMatrix] Window: %d SNPs around lead (idx %d, %d±%d of %d total)%n",
            windowBimIds.size(), leadIdx, leadIdx, ldWindow, allBimIds.size());

        // ── Extract window subset via PLINK --extract ──
        File extractFile = new File(ldDir, "_extract_window.txt");
        try (PrintWriter pw = new PrintWriter(new FileWriter(extractFile))) {
            for (String id : windowBimIds) pw.println(id);
        }

        String fullBfile = new File(matchedDir, "matched_ref").getAbsolutePath();
        String windowBfile = new File(ldDir, "_window_ref").getAbsolutePath();

        List<String> extractCmd = Arrays.asList(
            plinkBin, "--bfile", fullBfile,
            "--extract", extractFile.getAbsolutePath(),
            "--make-bed", "--out", windowBfile, "--silent"
        );
        try {
            Process proc = new ProcessBuilder(extractCmd).redirectErrorStream(true).start();
            try (BufferedReader drain = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                while (drain.readLine() != null) {}
            }
            proc.waitFor();
        } catch (InterruptedException e) { Thread.currentThread().interrupt(); }

        // ── Write SNP order for window ──
        try (PrintWriter pw = new PrintWriter(new FileWriter(new File(ldDir, "ld_snp_order.txt")))) {
            for (String snp : windowSnpKeys) pw.println(snp);
        }
        result.snpCount = windowSnpKeys.size();

        // ── Compute r matrix (signed correlation) ──
        boolean rOk = computeMatrix(plinkBin, windowBfile, ldDir, "r", "--r", "square");
        if (!rOk) {
            System.err.println("[LdMatrix] r matrix computation failed, trying r2 only");
        }

        // ── Compute r² matrix ──
        boolean r2Ok = computeMatrix(plinkBin, windowBfile, ldDir, "r2", "--r2", "square");
        if (!r2Ok) {
            result.error = "r² matrix computation failed";
            return result;
        }

        // ── Compute D' matrix ──
        boolean dpOk = computeDprime(plinkBin, windowBfile, ldDir, windowSnpKeys.size());
        if (!dpOk) {
            System.err.println("[LdMatrix] D' computation failed (optional)");
        }

        // Cleanup window temp files
        for (String ext : new String[]{".bed",".bim",".fam",".log",".nosex"}) {
            new File(windowBfile + ext).delete();
        }
        extractFile.delete();

        // Write manifest
        StepManifest manifest = new StepManifest();
        manifest.stepName = "ld";
        manifest.inputHash = inputHash;
        manifest.builtAt = System.currentTimeMillis();
        manifest.outputs.put("snp_count", String.valueOf(windowSnpKeys.size()));
        manifest.outputs.put("ld_window", String.valueOf(ldWindow));
        manifest.outputs.put("has_r", String.valueOf(rOk));
        manifest.outputs.put("has_r2", String.valueOf(r2Ok));
        manifest.outputs.put("has_dprime", String.valueOf(dpOk));
        manifest.write(ldDir);

        result.ok = true;
        System.out.printf("[LdMatrix] %d SNPs (window ±%d around lead): r=%s r2=%s dprime=%s%n",
            windowSnpKeys.size(), ldWindow, rOk, r2Ok, dpOk);
        return result;
    }

    private static boolean computeMatrix(String plinkBin, String bfilePrefix,
                                          File ldDir, String name,
                                          String plinkFlag, String shape) throws IOException {
        String tmpPrefix = new File(ldDir, "_plink_" + name).getAbsolutePath();

        List<String> cmd = Arrays.asList(
            plinkBin,
            "--bfile", bfilePrefix,
            plinkFlag, shape,
            "--out", tmpPrefix,
            "--silent"
        );

        try {
            Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                while (br.readLine() != null) {}
            }
            int exit = proc.waitFor();
            if (exit != 0) return false;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return false;
        }

        // PLINK writes .ld file — rename to our canonical name
        File plinkOut = new File(tmpPrefix + ".ld");
        File target = new File(ldDir, "ld_" + name + ".matrix");
        if (plinkOut.exists()) {
            target.delete();
            plinkOut.renameTo(target);
        }

        // Cleanup
        new File(tmpPrefix + ".log").delete();
        new File(tmpPrefix + ".nosex").delete();

        return target.exists();
    }

    private static final int DPRIME_MAX_SNPS = 5000;

    private static boolean computeDprime(String plinkBin, String bfilePrefix,
                                          File ldDir, int snpCount) throws IOException {
        // Skip D' for large loci — the n×n matrix doesn't fit in memory
        if (snpCount > DPRIME_MAX_SNPS) {
            System.out.printf("[LdMatrix] Skipping D' for %d SNPs (limit %d)%n",
                snpCount, DPRIME_MAX_SNPS);
            return false;
        }

        // Use PLINK --r2 dprime with a window, write sparse long-format output
        String tmpPrefix = new File(ldDir, "_plink_dprime").getAbsolutePath();
        int window = Math.min(snpCount, 500);

        List<String> cmd = Arrays.asList(
            plinkBin,
            "--bfile", bfilePrefix,
            "--r2", "dprime",
            "--ld-window", String.valueOf(window),
            "--ld-window-kb", "2000",
            "--ld-window-r2", "0",
            "--out", tmpPrefix,
            "--silent"
        );

        try {
            Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                while (br.readLine() != null) {}
            }
            int exit = proc.waitFor();
            if (exit != 0) return false;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return false;
        }

        File plinkOut = new File(tmpPrefix + ".ld");
        if (!plinkOut.exists()) return false;

        // Write D' as sparse long-format TSV (not dense matrix) to avoid OOM
        File target = new File(ldDir, "ld_dprime.matrix");
        try (BufferedReader br = new BufferedReader(new FileReader(plinkOut));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(target)))) {

            pw.println("SNP_A\tBP_A\tSNP_B\tBP_B\tR2\tDP");
            String line = br.readLine(); // skip PLINK header
            while ((line = br.readLine()) != null) {
                String[] f = line.trim().split("\\s+");
                if (f.length < 7) continue;
                String snpA = f[2], bpA = f[1];
                String snpB = f[5], bpB = f[4];
                String r2 = f[6];
                String dp = f.length > 7 ? f[7] : "NA";
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s%n", snpA, bpA, snpB, bpB, r2, dp);
            }
        }

        plinkOut.delete();
        new File(tmpPrefix + ".log").delete();
        new File(tmpPrefix + ".nosex").delete();

        return target.exists();
    }

    private static int findLeadSnpIndex(File harmonizedDir, List<Long> bimPositions) {
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists()) return -1;

        long bestPos = -1;
        double bestNlp = -1;
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            br.readLine(); // skip header
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;
                long pos = Long.parseLong(f[2]);
                double pval = Double.parseDouble(f[5]);
                if (pval <= 0 || Double.isNaN(pval)) continue;
                double nlp = -Math.log10(pval);
                if (nlp > bestNlp) { bestNlp = nlp; bestPos = pos; }
            }
        } catch (Exception e) { return -1; }

        if (bestPos < 0) return -1;
        for (int i = 0; i < bimPositions.size(); i++) {
            if (bimPositions.get(i) == bestPos) return i;
        }
        // Closest position if exact match not found
        int closest = 0;
        long closestDist = Long.MAX_VALUE;
        for (int i = 0; i < bimPositions.size(); i++) {
            long d = Math.abs(bimPositions.get(i) - bestPos);
            if (d < closestDist) { closestDist = d; closest = i; }
        }
        return closest;
    }

    private static double parseDouble(String s) {
        try {
            if (s == null || s.isEmpty() || s.equals("NA") || s.equals("nan")) return Double.NaN;
            return Double.parseDouble(s);
        } catch (NumberFormatException e) { return Double.NaN; }
    }
}
