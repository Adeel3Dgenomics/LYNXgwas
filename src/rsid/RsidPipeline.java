package rsid;

import java.io.*;
import java.nio.file.*;
import java.time.*;
import java.time.format.*;
import java.util.*;

/**
 * Orchestrates the full rsID recovery pipeline with progress tracking.
 * Uses the existing RsidRecovery/RsidMatcher module unchanged.
 */
public class RsidPipeline {

    /**
     * Run the full rsID recovery for a project.
     * @param projectDir  project directory (projects/{id})
     * @param gwasFile    path to GWAS file
     * @param lociFile    path to loci file
     * @param snpDbFolder path to dbSNP VCF folder
     * @param build       genome build (e.g. "hg19")
     * @param colChr      GWAS column name for chromosome
     * @param colPos      GWAS column name for position
     * @param colEa       GWAS column name for effect allele
     * @param colNea      GWAS column name for non-effect allele
     * @param progress    progress tracker (polled by UI)
     */
    public static void run(String projectDir, String gwasFile, String lociFile,
                           String snpDbFolder, String build,
                           String colChr, String colPos, String colEa, String colNea,
                           RsidProgress progress) {
        try {
            // Step 1: Validate inputs
            progress.stepIndex = 1;
            progress.currentStep = "Validating project & build";

            if (!new File(gwasFile).exists()) throw new IOException("GWAS file not found: " + gwasFile);
            if (!new File(lociFile).exists()) throw new IOException("Loci file not found: " + lociFile);
            if (!new File(snpDbFolder).exists()) throw new IOException("SNP database folder not found: " + snpDbFolder);

            // Parse loci
            List<String> lociChrs = new ArrayList<>();
            List<int[]> loci = new ArrayList<>();
            try (BufferedReader br = new BufferedReader(new FileReader(lociFile))) {
                String header = br.readLine();
                String[] hcols = header.trim().split("\t");
                int iChr = -1, iStart = -1, iEnd = -1;
                for (int i = 0; i < hcols.length; i++) {
                    String h = hcols[i].trim().toLowerCase();
                    if (h.equals("meta_chr"))   iChr = i;
                    if (h.equals("meta_start")) iStart = i;
                    if (h.equals("meta_end"))   iEnd = i;
                }
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.trim().isEmpty()) continue;
                    String[] f = line.trim().split("\t");
                    lociChrs.add(f[iChr].trim());
                    loci.add(new int[]{Integer.parseInt(f[iStart].trim()), Integer.parseInt(f[iEnd].trim())});
                }
            }

            progress.totalLoci = loci.size();
            System.out.printf("[RsidPipeline] %d loci, build=%s%n", loci.size(), build);

            // Verify needed chromosome VCFs exist
            Set<String> neededChrs = new LinkedHashSet<>(lociChrs);
            List<String> missingChrs = new ArrayList<>();
            for (String chr : neededChrs) {
                if (!new File(snpDbFolder, chr + ".vcf.gz").exists()) missingChrs.add(chr);
            }
            if (!missingChrs.isEmpty()) {
                System.err.printf("[RsidPipeline] Warning: missing VCFs for chromosomes: %s%n", missingChrs);
            }

            if (progress.cancelled) return;

            // Step 2: Build position indices from dbSNP
            progress.stepIndex = 2;
            progress.currentStep = "Loading dbSNP regions";

            RsidRecovery recovery = new RsidRecovery(snpDbFolder, build);
            Map<String, Map<Integer, List<DbSnpRecord>>> chrPosIndices = new LinkedHashMap<>();

            for (int li = 0; li < loci.size(); li++) {
                if (progress.cancelled) return;
                progress.currentLocus = li + 1;
                String chr = lociChrs.get(li);
                int start = loci.get(li)[0], end = loci.get(li)[1];
                String key = chr + ":" + start + "-" + end;

                Map<Integer, List<DbSnpRecord>> posIndex = recovery.buildPositionIndex(chr, start, end);
                chrPosIndices.put(key, posIndex);
                System.out.printf("[RsidPipeline] Locus %d/%d: %s — %,d dbSNP positions%n",
                    li + 1, loci.size(), key, posIndex.size());
            }

            if (progress.cancelled) return;

            // Steps 3-5: Stream GWAS file — match and write in one pass (no full-file buffering)
            progress.stepIndex = 3;
            progress.currentStep = "Matching & writing";
            progress.currentLocus = 0;

            // Always read from the ORIGINAL gwas file (strip .with_rsid if present)
            String originalGwas = gwasFile.replaceAll("\\.with_rsid\\.(tsv|txt|csv)$", "")
                                          .replaceAll("$", "");
            // If stripping made it lose the extension, it was already the original
            if (!new File(originalGwas).exists()) originalGwas = gwasFile;
            // Reconstruct: find the actual original with extension
            if (!new File(originalGwas).exists()) {
                for (String ext : new String[]{".tsv", ".txt", ".csv"}) {
                    String try_ = gwasFile.replaceAll("\\.with_rsid\\.(tsv|txt|csv)$", ext);
                    if (new File(try_).exists()) { originalGwas = try_; break; }
                }
            }
            if (!new File(originalGwas).exists()) {
                throw new IOException("Original GWAS file not found. Tried: " + originalGwas);
            }

            String baseName = originalGwas.replaceAll("\\.(tsv|txt|csv)$", "");
            String outputName = baseName + ".with_rsid.tsv";
            String csvPath = projectDir + "/rsid_recovery_results.csv";

            System.out.printf("[RsidPipeline] Reading original: %s%n", originalGwas);
            System.out.printf("[RsidPipeline] Writing to: %s%n", outputName);

            try (BufferedReader br = new BufferedReader(new FileReader(originalGwas), 1024 * 1024);
                 PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outputName), 1024 * 1024));
                 PrintWriter csv = new PrintWriter(new BufferedWriter(new FileWriter(csvPath)))) {

                String gwasHeader = br.readLine();
                if (gwasHeader == null) throw new IOException("GWAS file is empty");
                out.println(gwasHeader + "\trecovered_rsid");
                csv.println(MatchResult.csvHeader());

                String[] hcols = gwasHeader.trim().split("\t");
                int iChr = colIdx(hcols, colChr), iPos = colIdx(hcols, colPos);
                int iEa = colIdx(hcols, colEa), iNea = colIdx(hcols, colNea);
                long lineCount = 0;

                String line;
                while ((line = br.readLine()) != null) {
                    if (progress.cancelled) break;
                    if (line.isEmpty()) { out.println(line + "\t"); continue; }

                    lineCount++;
                    String chr = null;
                    int pos = -1;
                    MatchResult result = null;

                    int tabCount = 0;
                    for (int i = 0; i < line.length() && tabCount <= Math.max(Math.max(iChr, iPos), Math.max(iEa, iNea)); i++)
                        if (line.charAt(i) == '\t') tabCount++;

                    if (tabCount > Math.max(Math.max(iChr, iPos), Math.max(iEa, iNea))) {
                        String[] f = line.split("\t", -1);
                        chr = f[iChr].trim();
                        try { pos = Integer.parseInt(f[iPos].trim()); }
                        catch (NumberFormatException ignored) {}

                        if (pos >= 0) {
                            String ea = f[iEa].trim(), nea = f[iNea].trim();
                            for (int li = 0; li < loci.size(); li++) {
                                if (chr.equals(lociChrs.get(li))
                                        && pos >= loci.get(li)[0] && pos <= loci.get(li)[1]) {
                                    String key = chr + ":" + loci.get(li)[0] + "-" + loci.get(li)[1];
                                    Map<Integer, List<DbSnpRecord>> posIndex = chrPosIndices.get(key);
                                    if (posIndex != null)
                                        result = RsidMatcher.match(chr, pos, ea, nea, posIndex);
                                    break;
                                }
                            }
                        }
                    }

                    String rsid = (result != null && result.assignedRsid != null) ? result.assignedRsid : "";
                    out.println(line + "\t" + rsid);

                    if (result != null) {
                        csv.println(result.toCsv());
                        progress.totalSnps++;
                        if (result.assignedRsid != null) {
                            progress.matched.incrementAndGet();
                            if ("matched_forward".equals(result.matchReason))
                                progress.forward.incrementAndGet();
                            else
                                progress.reverse.incrementAndGet();
                        } else {
                            progress.unmatched.incrementAndGet();
                        }
                    }

                    if (lineCount % 2_000_000 == 0) {
                        System.out.printf("[RsidPipeline] %.0fM lines, %d matched%n",
                            lineCount / 1e6, progress.matched.get());
                    }
                }
            }

            if (progress.cancelled) return;

            progress.stepIndex = 4;
            progress.currentStep = "NCBI fallback (skipped)";

            // Compute final stats
            int totalMatched = progress.matched.get();
            int totalProcessed = progress.totalSnps;
            double rate = totalProcessed > 0 ? totalMatched * 100.0 / totalProcessed : 0;

            progress.stepIndex = 5;
            progress.currentStep = "Updating project";

            System.out.printf("%nTotal GWAS SNPs in loci: %d%n", totalProcessed);
            System.out.printf("Matched: %d (%.1f%%) — forward: %d, reverse: %d%n",
                totalMatched, rate, progress.forward.get(), progress.reverse.get());
            System.out.printf("Unmatched: %d%n", progress.unmatched.get());

            // ── Wire recovered rsIDs back into the project ──────────────
            // 1. Update config.properties: point to new GWAS file + set rsid column
            File configFile = new File(projectDir, "config.properties");
            if (configFile.exists()) {
                Properties props = new Properties();
                // Read with backslash normalization (same as Config.loadProperties)
                StringBuilder sb = new StringBuilder();
                try (BufferedReader cfgBr = new BufferedReader(new FileReader(configFile))) {
                    String cfgLine;
                    while ((cfgLine = cfgBr.readLine()) != null) {
                        if (!cfgLine.startsWith("#") && !cfgLine.startsWith("!") && cfgLine.contains("=")) {
                            int eq = cfgLine.indexOf('=');
                            String val = cfgLine.substring(eq + 1).replace("\\", "/");
                            cfgLine = cfgLine.substring(0, eq + 1) + val;
                        }
                        sb.append(cfgLine).append('\n');
                    }
                }
                props.load(new java.io.StringReader(sb.toString()));

                // Normalize the output path to forward slashes
                String normalizedOutput = outputName.replace("\\", "/");
                props.setProperty("gwas.file", normalizedOutput);
                props.setProperty("col.rsid", "recovered_rsid");

                // Rewrite config.properties with updated gwas.file and col.rsid
                // Read original lines, replace only the two keys we need to change
                List<String> configLines = Files.readAllLines(configFile.toPath());
                try (PrintWriter cpw = new PrintWriter(new BufferedWriter(new FileWriter(configFile)))) {
                    boolean wroteGwas = false, wroteRsid = false;
                    for (String cl : configLines) {
                        String trimCl = cl.trim();
                        if (trimCl.startsWith("gwas.file=")) {
                            cpw.println("gwas.file=" + normalizedOutput);
                            wroteGwas = true;
                        } else if (trimCl.startsWith("col.rsid=")) {
                            cpw.println("col.rsid=recovered_rsid");
                            wroteRsid = true;
                        } else {
                            cpw.println(cl);
                        }
                    }
                    if (!wroteGwas) cpw.println("gwas.file=" + normalizedOutput);
                    if (!wroteRsid) cpw.println("col.rsid=recovered_rsid");
                }

                System.out.printf("[RsidPipeline] Updated config.properties: gwas.file → %s, col.rsid → recovered_rsid%n", normalizedOutput);
            }

            // 2. Update project.json with rsID status
            try {
                File pjFile = new File(projectDir, "project.json");
                if (pjFile.exists()) {
                    String pjContent = new String(Files.readAllBytes(pjFile.toPath()), "UTF-8");
                    String dateStr = Instant.now().atOffset(ZoneOffset.UTC)
                        .format(DateTimeFormatter.ISO_OFFSET_DATE_TIME);
                    if (pjContent.contains("\"rsid_column_present\"")) {
                        pjContent = pjContent
                            .replaceAll("\"rsid_column_present\":\\s*\\w+", "\"rsid_column_present\": true")
                            .replaceAll("\"rsid_recovery_status\":\\s*\"[^\"]*\"", "\"rsid_recovery_status\": \"completed\"")
                            .replaceAll("\"rsid_recovery_rate\":\\s*\"[^\"]*\"", "\"rsid_recovery_rate\": \"" + String.format("%.1f%%", rate) + "\"")
                            .replaceAll("\"rsid_recovery_date\":\\s*\"[^\"]*\"", "\"rsid_recovery_date\": \"" + dateStr + "\"");
                    }
                    Files.writeString(pjFile.toPath(), pjContent);
                }
            } catch (Exception e) {
                System.err.println("[RsidPipeline] Could not update project.json: " + e.getMessage());
            }

            // 3. Clear fingerprint cache so the project is detected as stale
            //    (config.properties changed → core_input_fingerprint changed → reprocess triggers)
            Files.deleteIfExists(new File(projectDir, ".fingerprint_cache").toPath());

            progress.recoveryRate = rate;
            progress.outputFile = outputName;
            progress.done = true;
            progress.currentStep = "Complete";

            System.out.printf("[RsidPipeline] Done. Recovery rate: %.1f%%, output: %s%n", rate, outputName);
            System.out.println("[RsidPipeline] Project config updated — will reprocess on next run to load rsIDs into viewer.");

        } catch (Exception e) {
            progress.error = e.getMessage();
            progress.done = true;
            progress.currentStep = "Error";
            System.err.println("[RsidPipeline] Failed: " + e.getMessage());
            e.printStackTrace(System.err);
        }
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }
}
