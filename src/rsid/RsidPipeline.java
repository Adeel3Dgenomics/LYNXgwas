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
        run(projectDir, gwasFile, lociFile, snpDbFolder, build,
            colChr, colPos, colEa, colNea, progress, null);
    }

    public static void run(String projectDir, String gwasFile, String lociFile,
                           String snpDbFolder, String build,
                           String colChr, String colPos, String colEa, String colNea,
                           RsidProgress progress, RsidApiCompleter.Config apiConfig) {
        try {
            progress.totalSteps = 6;

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

            // ── Step 4: API completion for unmatched SNPs (optional) ──
            progress.stepIndex = 4;
            int localMatched = progress.matched.get();
            int localUnmatched = progress.unmatched.get();
            int totalProcessed = progress.totalSnps;

            System.out.printf("%n[RsidPipeline] Local matching done: %,d matched, %,d unmatched%n",
                localMatched, localUnmatched);

            Map<String, String> apiRecovered = new HashMap<>();
            if (apiConfig != null && apiConfig.enabled && localUnmatched > 0) {
                progress.currentStep = "API completion";
                // Collect unmatched SNPs from CSV
                List<RsidApiCompleter.UnmatchedSnp> unmatchedList = new ArrayList<>();
                try (BufferedReader csvBr2 = new BufferedReader(new FileReader(csvPath))) {
                    csvBr2.readLine(); // skip header
                    String csvLine2;
                    while ((csvLine2 = csvBr2.readLine()) != null) {
                        String[] parts = csvLine2.split(",", 7);
                        if (parts.length >= 6 && (parts[5].equals("no_pos") || parts[5].equals("pos_only_allele_mismatch"))) {
                            RsidApiCompleter.UnmatchedSnp us = new RsidApiCompleter.UnmatchedSnp();
                            us.chr = parts[0]; us.pos = Long.parseLong(parts[1]);
                            us.ea = parts[2]; us.nea = parts[3];
                            unmatchedList.add(us);
                        }
                    }
                }

                RsidApiCache cache = new RsidApiCache(projectDir);
                RsidApiCompleter.CompletionResult apiResult =
                    RsidApiCompleter.complete(unmatchedList, build, apiConfig, cache, progress);
                apiRecovered = apiResult.recovered;

                progress.matched.addAndGet(apiResult.apiMatched);
                progress.unmatched.addAndGet(-apiResult.apiMatched);
                System.out.printf("[RsidPipeline] Combined recovery: %.2f%%  (%,d / %,d)%n",
                    (localMatched + apiResult.apiMatched) * 100.0 / totalProcessed,
                    localMatched + apiResult.apiMatched, totalProcessed);
            } else if (apiConfig != null && !apiConfig.enabled) {
                progress.currentStep = "API completion (disabled)";
                System.out.println("[RsidPipeline] API completion: disabled by user");
            } else {
                progress.currentStep = "API completion (skipped)";
            }

            if (progress.cancelled) return;

            // ── Step 5: Cross-file GWAS lookup (Strategy 3) ──
            // Checks every project's GWAS file that has a genuine rsID column (including
            // this project's own file) for the same chr:pos, in case another cohort's
            // summary stats already resolved it. Runs after local dbSNP + API completion;
            // whatever's still unmatched afterward goes to the cross-project
            // "missing rsIDs" list (see LocalServer /api/missing-rsids).
            progress.stepIndex = 5;
            Map<String, String> crossFileRecovered = new HashMap<>();
            if (progress.unmatched.get() > 0) {
                progress.currentStep = "Cross-file GWAS lookup";
                Map<String, String> crossIndex = CrossFileLookup.buildIndex(new File(projectDir).getParentFile());
                try (BufferedReader csvBr3 = new BufferedReader(new FileReader(csvPath))) {
                    csvBr3.readLine(); // skip header
                    String csvLine3;
                    while ((csvLine3 = csvBr3.readLine()) != null) {
                        String[] parts = csvLine3.split(",", 7);
                        if (parts.length < 6) continue;
                        String key = parts[0] + ":" + parts[1];
                        boolean alreadyResolved = !parts[4].isEmpty() || apiRecovered.containsKey(key);
                        if (alreadyResolved) continue;
                        String rsid = crossIndex.get(key);
                        if (rsid != null) crossFileRecovered.put(key, rsid);
                    }
                }
                progress.matched.addAndGet(crossFileRecovered.size());
                progress.unmatched.addAndGet(-crossFileRecovered.size());
                System.out.printf("[RsidPipeline] Cross-file lookup: %,d additional matches (index size %,d)%n",
                    crossFileRecovered.size(), crossIndex.size());
            } else {
                progress.currentStep = "Cross-file GWAS lookup (skipped — nothing unmatched)";
            }

            int totalMatched = progress.matched.get();
            double rate = totalProcessed > 0 ? totalMatched * 100.0 / totalProcessed : 0;

            progress.stepIndex = 6;
            progress.currentStep = "Patching locus JSONs";

            System.out.printf("Total matched: %,d (%.1f%%) — forward: %,d, reverse: %,d%n",
                totalMatched, rate, progress.forward.get(), progress.reverse.get());
            System.out.printf("Still unmatched: %,d%n", progress.unmatched.get());

            // ── Persist API + cross-file recoveries back into the CSV so the cross-project
            //    "missing rsIDs" list stays accurate (previously only local matches were
            //    ever written back, so API-recovered SNPs would wrongly keep showing as missing) ──
            if (!apiRecovered.isEmpty() || !crossFileRecovered.isEmpty()) {
                List<String> csvLines = Files.readAllLines(Paths.get(csvPath));
                List<String> updated = new ArrayList<>(csvLines.size());
                if (!csvLines.isEmpty()) updated.add(csvLines.get(0)); // header
                for (int i = 1; i < csvLines.size(); i++) {
                    String[] parts = csvLines.get(i).split(",", -1);
                    if (parts.length >= 6 && parts[4].isEmpty()) {
                        String key = parts[0] + ":" + parts[1];
                        String rsid = apiRecovered.get(key);
                        String reason = "api";
                        if (rsid == null) { rsid = crossFileRecovered.get(key); reason = "cross_file"; }
                        if (rsid != null) { parts[4] = rsid; parts[5] = reason; }
                    }
                    updated.add(String.join(",", parts));
                }
                Files.write(Paths.get(csvPath), updated);
            }

            // ── Patch existing locus JSONs with recovered rsIDs (no reprocess needed) ──
            // Build chr:pos → rsid lookup from the (now fully updated) CSV
            Map<String, String> posToRsid = new HashMap<>();
            try (BufferedReader csvBr = new BufferedReader(new FileReader(csvPath))) {
                csvBr.readLine(); // skip header
                String csvLine;
                while ((csvLine = csvBr.readLine()) != null) {
                    String[] parts = csvLine.split(",", 6);
                    if (parts.length >= 5 && !parts[4].isEmpty()) {
                        posToRsid.put(parts[0] + ":" + parts[1], parts[4]); // chr:pos → rsid
                    }
                }
            }
            System.out.printf("[RsidPipeline] Patching locus JSONs with %,d rsIDs (local + API + cross-file)%n",
                posToRsid.size());

            int patchedCount = patchLocusJsonsWithRsids(projectDir, posToRsid);
            System.out.printf("[RsidPipeline] Patched %d locus files%n", patchedCount);

            // Update project.json with rsID status (don't touch config.properties — avoids stale fingerprint)
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

            progress.recoveryRate = rate;
            progress.outputFile = outputName;
            progress.done = true;
            progress.currentStep = "Complete";

            System.out.printf("[RsidPipeline] Done. Recovery rate: %.1f%%, output: %s%n", rate, outputName);
            System.out.println("[RsidPipeline] Done — rsIDs patched into locus JSONs. No reprocess needed.");

        } catch (Exception e) {
            progress.error = e.getMessage();
            progress.done = true;
            progress.currentStep = "Error";
            System.err.println("[RsidPipeline] Failed: " + e.getMessage());
            e.printStackTrace(System.err);
        }
    }

    /**
     * Patches "id" fields in a project's locus_N.json files wherever chr:pos matches
     * a key in posToRsid. Shared by the full recovery pipeline and by manual
     * single-SNP rsID submission (see LocalServer's /api/submit-rsid).
     * Returns the number of locus files actually changed.
     */
    public static int patchLocusJsonsWithRsids(String projectDir, Map<String, String> posToRsid) throws IOException {
        File dataDir = new File(projectDir, "data");
        if (!dataDir.isDirectory()) return 0;
        File[] locusFiles = dataDir.listFiles((d, n) -> n.matches("locus_\\d+\\.json"));
        if (locusFiles == null) return 0;

        int patched = 0;
        for (File lf : locusFiles) {
            String content = new String(Files.readAllBytes(lf.toPath()), "UTF-8");
            boolean changed = false;

            // Replace SNP IDs: find "id":"oldId","chr":"X","pos":NNN patterns
            // and replace oldId with rsid from our lookup
            StringBuilder result = new StringBuilder(content.length());
            int idx = 0;
            while (idx < content.length()) {
                int idStart = content.indexOf("\"id\":\"", idx);
                if (idStart < 0) { result.append(content, idx, content.length()); break; }
                result.append(content, idx, idStart);

                // Extract current id
                int valStart = idStart + 6;
                int valEnd = content.indexOf('"', valStart);
                if (valEnd < 0) { result.append(content, idStart, content.length()); break; }
                String oldId = content.substring(valStart, valEnd);

                // Look ahead for chr and pos
                int chrIdx = content.indexOf("\"chr\":\"", valEnd);
                int posIdx = content.indexOf("\"pos\":", valEnd);
                if (chrIdx >= 0 && posIdx >= 0 && chrIdx - valEnd < 100 && posIdx - valEnd < 150) {
                    int chrValStart = chrIdx + 7;
                    int chrValEnd = content.indexOf('"', chrValStart);
                    int posValStart = posIdx + 6;
                    int posValEnd = posValStart;
                    while (posValEnd < content.length() && Character.isDigit(content.charAt(posValEnd))) posValEnd++;

                    if (chrValEnd > chrValStart && posValEnd > posValStart) {
                        String chr = content.substring(chrValStart, chrValEnd);
                        String pos = content.substring(posValStart, posValEnd);
                        String rsid = posToRsid.get(chr + ":" + pos);
                        if (rsid != null && !rsid.equals(oldId)) {
                            result.append("\"id\":\"").append(rsid).append('"');
                            idx = valEnd + 1;
                            changed = true;
                            continue;
                        }
                    }
                }

                result.append(content, idStart, valEnd + 1);
                idx = valEnd + 1;
            }

            if (changed) {
                Files.writeString(lf.toPath(), result.toString());
                patched++;
            }
        }
        return patched;
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }
}
