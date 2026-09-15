
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * SuSiEx (Yuan et al.) cross-ancestry fine-mapping adapter.
 *
 * Unlike the single-project tools in tools/*.yaml (SuSiE, FINEMAP, COJO, coloc),
 * SuSiEx takes one summary-statistics + LD input per ancestry/cohort and jointly
 * fine-maps a shared locus, so it does not fit PluginEngine's one-project-one-locus
 * RunRequest model. It is driven instead by a dedicated multi-project endpoint
 * (see LocalServer#susiexRun) that gathers already-built base-pipeline artifacts
 * (harmonized_gwas.tsv + the locus's own matched_ref PLINK subset) from each
 * selected project/locus and hands them here.
 *
 * SuSiEx itself is an external Python tool (github.com/getian107/SuSiEx) that must
 * already be installed; this adapter only prepares its per-ancestry input files in
 * the exact column layout it is told to expect (via explicit --*_col indices, so
 * behavior does not depend on the installed version's default column-name guesses)
 * and invokes it, pointing --ld_file directly at each project's already-computed
 * matched_ref PLINK bfile prefix (SuSiEx computes LD itself via plink).
 *
 * Output parsing is intentionally schema-tolerant: SuSiEx's exact output column
 * names have varied across versions. We look for the "*.cs" (credible set) or
 * "*.snp"/"*.summary" file it writes, then locate a PIP-like and an id-like column
 * by name rather than a fixed position. If nothing recognizable is found, the raw
 * output directory listing and run.log are surfaced so the run isn't a silent
 * black box.
 */
public class SusiexAdapter {

    public static class Member {
        public String label;        // e.g. project name / ancestry
        public File harmonizedDir;  // .../harmonized/harmonized_gwas.tsv
        public File matchedDir;     // .../matched/matched_ref.{bed,bim,fam}
        public int sampleN;
    }

    public static class ResultRow {
        public String snpId;
        public String chr;
        public long pos = -1;
        public double pip = Double.NaN;
        public String csId = "";
    }

    public static class RunResult {
        public boolean ok;
        public String error;
        public List<ResultRow> rows = new ArrayList<>();
        public int nCredibleSets;
        public String rawOutputFile;
        public String logTail;
    }

    public static RunResult prepareAndRun(List<Member> members, File runDir,
                                           String chr, long start, long end,
                                           String susiexPath, String plinkPath,
                                           int maxCausal, double pvalThresh, double mafThresh) throws IOException {
        if (members.size() < 2)
            throw new IOException("SuSiEx needs at least two ancestry/cohort inputs.");
        runDir.mkdirs();

        List<String> sstFiles = new ArrayList<>();
        List<String> ldFiles = new ArrayList<>();
        List<String> nGwas = new ArrayList<>();

        for (int i = 0; i < members.size(); i++) {
            Member m = members.get(i);
            if (m.sampleN <= 0)
                throw new IOException("Sample size (N) is required for member '" + m.label + "'.");
            File bim = new File(m.matchedDir, "matched_ref.bim");
            if (!bim.exists())
                throw new IOException("matched_ref.bim not found for '" + m.label + "'. Run the base pipeline for that locus first.");

            File sst = new File(runDir, "sumstats_" + (i + 1) + ".tsv");
            int n = writeSumstats(m.harmonizedDir, sst);
            if (n == 0)
                throw new IOException("No usable SNPs (beta/se) for member '" + m.label + "'.");

            sstFiles.add(sst.getAbsolutePath().replace("\\", "/"));
            ldFiles.add(new File(m.matchedDir, "matched_ref").getAbsolutePath().replace("\\", "/"));
            nGwas.add(String.valueOf(m.sampleN));
        }

        String outName = "result";
        List<String> cmd = new ArrayList<>();
        cmd.add("python");
        cmd.add(susiexPath);
        cmd.add("--sst_file=" + String.join(",", sstFiles));
        cmd.add("--n_gwas=" + String.join(",", nGwas));
        cmd.add("--ld_file=" + String.join(",", ldFiles));
        cmd.add("--out_dir=" + runDir.getAbsolutePath());
        cmd.add("--out_name=" + outName);
        cmd.add("--chr=" + chr);
        cmd.add("--bp=" + start + "," + end);
        cmd.add("--plink=" + (plinkPath != null && !plinkPath.isEmpty() ? plinkPath : "plink"));
        cmd.add("--rsid_col=1");
        cmd.add("--chr_col=2");
        cmd.add("--bp_col=3");
        cmd.add("--a1_col=4");
        cmd.add("--a2_col=5");
        cmd.add("--eff_col=6");
        cmd.add("--se_col=7");
        cmd.add("--pval_col=8");
        cmd.add("--pval_thresh=" + pvalThresh);
        cmd.add("--maf_thresh=" + mafThresh);
        cmd.add("--n_sig=" + maxCausal);
        cmd.add("--level=0.95");

        File logFile = new File(runDir, "run.log");
        System.out.printf("[SusiexAdapter] Running: %s%n", String.join(" ", cmd));

        RunResult result = new RunResult();
        int exitCode;
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd).directory(runDir).redirectErrorStream(true);
            Process proc = pb.start();
            List<String> tail = new ArrayList<>();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                 PrintWriter log = new PrintWriter(new BufferedWriter(new FileWriter(logFile)))) {
                String line;
                while ((line = br.readLine()) != null) {
                    log.println(line);
                    tail.add(line);
                    if (tail.size() > 60) tail.remove(0);
                }
            }
            boolean finished = proc.waitFor(30, TimeUnit.MINUTES);
            if (!finished) { proc.destroyForcibly(); throw new IOException("SuSiEx timed out after 30 minutes"); }
            exitCode = proc.exitValue();
            result.logTail = String.join("\n", tail);
        } catch (IOException e) {
            throw e;
        } catch (Exception e) {
            throw new IOException("Failed to launch SuSiEx: " + e.getMessage(), e);
        }

        if (exitCode != 0) {
            result.error = "SuSiEx exited with code " + exitCode + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }

        parseOutput(runDir, outName, result);
        result.ok = !result.rows.isEmpty();
        if (!result.ok && result.error == null) {
            result.error = "SuSiEx ran successfully but no recognizable output was found in "
                + runDir.getAbsolutePath() + ". Files present: " + Arrays.toString(runDir.list());
        }
        return result;
    }

    /** rsid  chr  bp  a1  a2  beta  se  pval (fixed column order, matches the --*_col indices passed to SuSiEx). */
    private static int writeSumstats(File harmonizedDir, File out) throws IOException {
        File gwas = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwas.exists()) throw new IOException("harmonized_gwas.tsv not found in " + harmonizedDir);
        int count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(gwas));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)))) {
            pw.println("rsid\tchr\tbp\ta1\ta2\tbeta\tse\tpval");
            br.readLine(); // snp_id chr pos ea nea pvalue beta se or n maf info rsid varid
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 8) continue;
                String beta = f[6], se = f[7], pval = f[5];
                if (beta.equals("NA") || se.equals("NA")) continue;
                String id = (f.length > 12 && !f[12].equals("NA") && !f[12].trim().isEmpty()) ? f[12] : f[0];
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n", id, f[1], f[2], f[3], f[4], beta, se, pval);
                count++;
            }
        }
        return count;
    }

    private static void parseOutput(File runDir, String outName, RunResult result) {
        File[] candidates = {
            new File(runDir, outName + ".cs"),
            new File(runDir, outName + ".snp"),
            new File(runDir, outName + ".summary"),
        };
        File found = null;
        for (File f : candidates) if (f.exists() && f.length() > 0) { found = f; break; }
        if (found == null) {
            // Fall back to any file this run just produced that isn't an input we wrote ourselves.
            File[] all = runDir.listFiles((d, n) -> n.startsWith(outName) && !n.equals("run.log"));
            if (all != null && all.length > 0) found = all[0];
        }
        if (found == null) return;
        result.rawOutputFile = found.getAbsolutePath();

        try (BufferedReader br = new BufferedReader(new FileReader(found))) {
            String header = br.readLine();
            if (header == null) return;
            String delim = header.contains("\t") ? "\t" : "\\s+";
            String[] cols = header.split(delim, -1);

            int idIdx = -1, pipIdx = -1, csIdx = -1, chrIdx = -1, posIdx = -1;
            for (int i = 0; i < cols.length; i++) {
                String c = cols[i].trim().toLowerCase();
                if (idIdx < 0 && (c.equals("snp") || c.contains("rsid") || c.contains("variant"))) idIdx = i;
                if (pipIdx < 0 && c.contains("pip")) pipIdx = i;
                if (csIdx < 0 && (c.equals("cs") || c.contains("cs_id") || c.equals("credible_set"))) csIdx = i;
                if (chrIdx < 0 && c.equals("chr")) chrIdx = i;
                if (posIdx < 0 && (c.equals("bp") || c.equals("pos"))) posIdx = i;
            }
            if (idIdx < 0 || pipIdx < 0) return; // unrecognized schema; caller reports the raw file path

            Set<String> csIds = new LinkedHashSet<>();
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] f = line.split(delim, -1);
                if (f.length <= Math.max(idIdx, pipIdx)) continue;
                ResultRow row = new ResultRow();
                row.snpId = f[idIdx].trim();
                try { row.pip = Double.parseDouble(f[pipIdx].trim()); } catch (NumberFormatException e) { continue; }
                if (chrIdx >= 0 && f.length > chrIdx) row.chr = f[chrIdx].trim();
                if (posIdx >= 0 && f.length > posIdx) {
                    try { row.pos = Long.parseLong(f[posIdx].trim()); } catch (NumberFormatException ignored) {}
                }
                if (csIdx >= 0 && f.length > csIdx) { row.csId = f[csIdx].trim(); csIds.add(row.csId); }
                result.rows.add(row);
            }
            result.nCredibleSets = csIds.size();
        } catch (IOException e) {
            result.error = "Failed to parse SuSiEx output " + found.getAbsolutePath() + ": " + e.getMessage();
        }
    }
}
