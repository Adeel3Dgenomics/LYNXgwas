package loci;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Identifies GWAS loci using the correct pipeline:
 * 1. Stream GWAS, pre-filter to p <= preFilterP (5e-3), group by chromosome
 * 2. Ref-panel match by chr:pos — drop variants absent from ref panel BIM
 * 3. Seed filter: p <= leadPThreshold
 * 4. PLINK clumping per chromosome (--clump-r2 0.6, --clump-kb 250)
 * 5. Merge clump intervals within 250kb into genomic loci
 * 6. Assign global indices, write output
 */
public class LociIdentifier {

    public static class Params {
        public double preFilterP       = 5e-3;
        public double leadPThreshold   = 5e-8;
        public double clumpR2          = 0.6;
        public int    clumpKb          = 250;
        public int    mergeDistanceBp  = 250_000;
        public int    minSnpsPerLocus  = 1;
        public String refPanelPath     = "";   // PLINK bfile prefix (e.g. /path/to/g1000_eas)
        public String plinkBin         = "";   // path to plink binary
        public int    parallelJobs     = 4;
    }

    public static class IdentifiedLocus {
        public int    globalIndex;
        public String chr;
        public long   start;
        public long   end;
        public String topSnpId;
        public double topP;
        public int    nSnps;

        public String toTsv() {
            return String.join("\t",
                String.valueOf(globalIndex), chr,
                String.valueOf(start), String.valueOf(end),
                String.valueOf(end - start + 1),
                topSnpId, String.format("%.4e", topP),
                String.valueOf(nSnps));
        }

        public static String tsvHeader() {
            return "genomic_locus\tmeta_chr\tmeta_start\tmeta_end\tsize_bp\ttop_snp\ttop_p\tn_snps";
        }
    }

    static class CandidateSnp {
        String id;       // rsID from ref panel (matched by position)
        String chr;
        long pos;
        double p;
    }

    // Clump interval parsed from PLINK .clumped output
    static class ClumpInterval {
        String chr;
        long start, end;
        String leadSnp;
        double leadP;
        int memberCount;
    }

    public static List<IdentifiedLocus> identify(
            String gwasFile, String colChr, String colPos, String colPvalue,
            String colId, Params params, LociProgress progress) throws IOException {

        // ── Stage A: Stream GWAS, pre-filter p <= preFilterP, group by chr ──
        if (progress != null) { progress.stepIndex = 1; progress.currentStep = "Pre-filtering GWAS (p <= " + params.preFilterP + ")"; }

        Map<String, List<CandidateSnp>> byChr = new TreeMap<>(LociIdentifier::chrCompare);
        long totalLines = 0, stageA = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) throw new IOException("GWAS file is empty");
            String[] cols = header.trim().split("\t");
            int iChr = colIdx(cols, colChr), iPos = colIdx(cols, colPos);
            int iP = colIdx(cols, colPvalue), iId = colIdx(cols, colId);
            if (iChr < 0 || iPos < 0 || iP < 0)
                throw new IOException("Missing required columns: chr=" + colChr + " pos=" + colPos + " p=" + colPvalue);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                totalLines++;
                String[] f = line.split("\t", -1);
                if (f.length <= Math.max(iChr, Math.max(iPos, iP))) continue;

                double p; long pos;
                try {
                    p = Double.parseDouble(f[iP].trim());
                    pos = Long.parseLong(f[iPos].trim());
                } catch (NumberFormatException e) { continue; }

                if (p <= 0 || p > params.preFilterP) continue;
                stageA++;

                String chr = f[iChr].trim();
                String id = (iId >= 0 && iId < f.length) ? f[iId].trim() : chr + ":" + pos;

                CandidateSnp snp = new CandidateSnp();
                snp.id = id; snp.chr = chr; snp.pos = pos; snp.p = p;
                byChr.computeIfAbsent(chr, k -> new ArrayList<>()).add(snp);

                if (totalLines % 5_000_000 == 0) {
                    System.out.printf("[LociIdentifier] %.0fM lines, %,d candidates%n", totalLines / 1e6, stageA);
                    if (progress != null) progress.candidateSnps = (int) stageA;
                }
            }
        }
        System.out.printf("[Loci] Stage A — pre-filter p<=%.0e: %,d / %,d lines%n", params.preFilterP, stageA, totalLines);
        if (progress != null) { progress.candidateSnps = (int) stageA; progress.totalChromosomes = byChr.size(); }

        // ── Stage B: Ref-panel match by chr:pos (per-chromosome BIM loading) ──
        if (progress != null) { progress.stepIndex = 2; progress.currentStep = "Ref-panel match (chr:pos)"; }
        long stageB = 0;

        if (!params.refPanelPath.isEmpty()) {
            for (Map.Entry<String, List<CandidateSnp>> entry : byChr.entrySet()) {
                String chr = entry.getKey();
                Set<Long> refPositions = loadBimPositions(params.refPanelPath, chr);
                if (refPositions.isEmpty()) {
                    System.out.printf("[Loci] chr%s: no BIM positions found, keeping all%n", chr);
                    stageB += entry.getValue().size();
                    continue;
                }
                int before = entry.getValue().size();

                // Match by chr:pos and replace ID with ref panel rsID
                Map<Long, String> posToRsid = loadBimPosToRsid(params.refPanelPath, chr);
                Iterator<CandidateSnp> it = entry.getValue().iterator();
                while (it.hasNext()) {
                    CandidateSnp snp = it.next();
                    if (!refPositions.contains(snp.pos)) {
                        it.remove();
                    } else {
                        String rsid = posToRsid.get(snp.pos);
                        if (rsid != null) snp.id = rsid;
                    }
                }
                int after = entry.getValue().size();
                stageB += after;
                System.out.printf("[Loci] chr%s: ref-panel match %,d → %,d%n", chr, before, after);
            }
            byChr.entrySet().removeIf(e -> e.getValue().isEmpty());
        } else {
            stageB = stageA;
            System.out.println("[Loci] Stage B — no ref panel, skipping match");
        }
        System.out.printf("[Loci] Stage B — after ref-panel match: %,d SNPs%n", stageB);

        // ── Stage C: Seed filter p <= leadPThreshold ──
        // Non-destructive: byChr must keep all preFilterP-level candidates so PLINK's
        // --clump-p2 can pull in suggestive (non-genome-wide-significant) SNPs as LD-linked
        // clump members. Only seedsByChr (used for counting and the no-PLINK fallback) is
        // restricted to lead-threshold-significant SNPs.
        if (progress != null) { progress.stepIndex = 3; progress.currentStep = "Seed filter (p <= " + params.leadPThreshold + ")"; }
        Map<String, List<CandidateSnp>> seedsByChr = new TreeMap<>(LociIdentifier::chrCompare);
        long stageC = 0;
        for (Map.Entry<String, List<CandidateSnp>> entry : byChr.entrySet()) {
            List<CandidateSnp> seeds = new ArrayList<>();
            for (CandidateSnp s : entry.getValue()) if (s.p <= params.leadPThreshold) seeds.add(s);
            if (!seeds.isEmpty()) { seedsByChr.put(entry.getKey(), seeds); stageC += seeds.size(); }
        }
        System.out.printf("[Loci] Stage C — seeds p<=%.0e: %,d across %d chromosomes%n",
            params.leadPThreshold, stageC, seedsByChr.size());
        if (progress != null) { progress.seedSnps = (int) stageC; }

        // Chromosomes with no seed SNP can't produce a clump at all — drop them from the
        // clumping input, but keep the full preFilterP-level pool for the ones that remain.
        byChr.keySet().retainAll(seedsByChr.keySet());

        // ── Stage D: PLINK clumping per chromosome ──
        List<ClumpInterval> allClumps;
        String plinkBin = params.plinkBin.isEmpty() ? findPlink() : params.plinkBin;

        if (plinkBin != null && !params.refPanelPath.isEmpty()) {
            if (progress != null) { progress.stepIndex = 4; progress.currentStep = "PLINK clumping"; }
            allClumps = runClumping(byChr, params, plinkBin, progress);
            System.out.printf("[Loci] Stage D — after PLINK clumping: %,d independent signals%n", allClumps.size());
        } else {
            // Fallback: no PLINK available — use positional merge directly on seeds.
            // No LD information is available here, so members can't be identified beyond
            // the seed SNP itself (unlike the PLINK path, which can pull in suggestive
            // LD-linked neighbors as members).
            if (plinkBin == null) System.out.println("[Loci] Stage D — PLINK not found, falling back to positional merge");
            else System.out.println("[Loci] Stage D — no ref panel for clumping, using positional merge");
            if (progress != null) { progress.stepIndex = 4; progress.currentStep = "Positional merge (no PLINK)"; }
            allClumps = new ArrayList<>();
            for (Map.Entry<String, List<CandidateSnp>> entry : seedsByChr.entrySet()) {
                String chr = entry.getKey();
                List<CandidateSnp> snps = entry.getValue();
                snps.sort(Comparator.comparingLong((CandidateSnp s) -> s.pos));
                for (CandidateSnp snp : snps) {
                    ClumpInterval ci = new ClumpInterval();
                    ci.chr = chr; ci.start = snp.pos; ci.end = snp.pos;
                    ci.leadSnp = snp.id; ci.leadP = snp.p; ci.memberCount = 1;
                    allClumps.add(ci);
                }
            }
        }

        // ── Stage E: Merge clump intervals within mergeDistanceBp ──
        if (progress != null) { progress.stepIndex = 5; progress.currentStep = "Merging clump intervals"; }
        List<IdentifiedLocus> allLoci = mergeClumpIntervals(allClumps, params.mergeDistanceBp, params.minSnpsPerLocus);
        System.out.printf("[Loci] Stage E — after 250kb clump merge: %,d loci%n", allLoci.size());

        // ── Stage F: Assign global indices ──
        if (progress != null) { progress.stepIndex = 6; progress.currentStep = "Finalizing"; }
        allLoci.sort(Comparator.<IdentifiedLocus>comparingInt(l -> chrToInt(l.chr)).thenComparingLong(l -> l.start));
        for (int i = 0; i < allLoci.size(); i++) allLoci.get(i).globalIndex = i + 1;

        System.out.printf("[Loci] Stage F — final: %d loci%n", allLoci.size());
        if (progress != null) {
            progress.lociFound = allLoci.size();
            progress.done = true;
            progress.currentStep = "Complete";
        }
        return allLoci;
    }

    // ── Per-chromosome BIM loading (by position, not by ID) ────────────

    static Set<Long> loadBimPositions(String plinkPrefix, String chr) {
        Set<Long> positions = new HashSet<>();
        // Try per-chromosome BIM first, then whole-genome BIM
        for (String suffix : new String[]{"_chr" + chr + ".bim", ".bim"}) {
            File bimFile = new File(plinkPrefix + suffix);
            if (!bimFile.exists()) continue;
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile), 256 * 1024)) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t");
                    if (f.length < 4) { f = line.split("\\s+"); }
                    if (f.length < 4) continue;
                    // BIM: chr, rsid, cm, pos, a1, a2
                    if (!f[0].equals(chr)) continue;
                    try { positions.add(Long.parseLong(f[3].trim())); }
                    catch (NumberFormatException ignored) {}
                }
            } catch (IOException e) {
                System.err.printf("[Loci] Error reading BIM for chr%s: %s%n", chr, e.getMessage());
            }
            break; // use first found
        }
        return positions;
    }

    static Map<Long, String> loadBimPosToRsid(String plinkPrefix, String chr) {
        Map<Long, String> map = new HashMap<>();
        for (String suffix : new String[]{"_chr" + chr + ".bim", ".bim"}) {
            File bimFile = new File(plinkPrefix + suffix);
            if (!bimFile.exists()) continue;
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile), 256 * 1024)) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t");
                    if (f.length < 4) { f = line.split("\\s+"); }
                    if (f.length < 4) continue;
                    if (!f[0].equals(chr)) continue;
                    try {
                        long pos = Long.parseLong(f[3].trim());
                        map.put(pos, f[1].trim());
                    } catch (NumberFormatException ignored) {}
                }
            } catch (IOException e) { break; }
            break;
        }
        return map;
    }

    // ── PLINK clumping ─────────────────────────────────────────────────

    static List<ClumpInterval> runClumping(Map<String, List<CandidateSnp>> byChr,
                                           Params params, String plinkBin,
                                           LociProgress progress) throws IOException {
        Path tmpDir = Files.createTempDirectory("lynx_clump_");
        List<ClumpInterval> allClumps = Collections.synchronizedList(new ArrayList<>());
        List<String> chrList = new ArrayList<>(byChr.keySet());
        int chromDone = 0;

        // Process chromosomes (sequentially to avoid heavy memory use)
        for (String chr : chrList) {
            if (progress != null && progress.cancelled) break;
            List<CandidateSnp> snps = byChr.get(chr);
            if (snps.isEmpty()) continue;

            // Write association file for PLINK --clump
            Path assocFile = tmpDir.resolve("chr" + chr + ".assoc");
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(assocFile.toFile())))) {
                pw.println("SNP\tP");
                for (CandidateSnp s : snps) pw.printf("%s\t%.6e%n", s.id, s.p);
            }

            // Run PLINK clump
            Path outPrefix = tmpDir.resolve("clump_chr" + chr);
            List<String> cmd = Arrays.asList(
                plinkBin,
                "--bfile", params.refPanelPath,
                "--clump", assocFile.toString(),
                "--clump-snp-field", "SNP",
                "--clump-field", "P",
                "--clump-p1", String.valueOf(params.leadPThreshold),
                "--clump-p2", String.valueOf(params.preFilterP),
                "--clump-r2", String.valueOf(params.clumpR2),
                "--clump-kb", String.valueOf(params.clumpKb),
                "--chr", chr,
                "--out", outPrefix.toString(),
                "--silent"
            );

            System.out.printf("[Loci] PLINK clump chr%s: %d seeds%n", chr, snps.size());
            try {
                Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
                // Drain output
                try (BufferedReader pbr = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                    while (pbr.readLine() != null) {}
                }
                int exit = proc.waitFor();
                if (exit != 0) {
                    System.err.printf("[Loci] PLINK clump chr%s exit code %d%n", chr, exit);
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }

            // Parse .clumped output
            Path clumpedFile = Paths.get(outPrefix + ".clumped");
            if (Files.exists(clumpedFile)) {
                List<ClumpInterval> chrClumps = parseClumpedFile(clumpedFile, chr, snps);
                allClumps.addAll(chrClumps);
                System.out.printf("[Loci] chr%s: %d clumps%n", chr, chrClumps.size());
            } else {
                System.out.printf("[Loci] chr%s: no .clumped output (no clumps)%n", chr);
            }

            chromDone++;
            if (progress != null) {
                progress.currentChromosome = chromDone;
                progress.lociFound = allClumps.size();
            }
        }

        // Cleanup temp files
        try {
            Files.walk(tmpDir).sorted(Comparator.reverseOrder()).map(Path::toFile).forEach(File::delete);
        } catch (IOException ignored) {}

        return allClumps;
    }

    /**
     * Parse PLINK .clumped output. Each row is an independent signal.
     * Build clump intervals from lead SNP position ± member positions.
     */
    static List<ClumpInterval> parseClumpedFile(Path clumpedFile, String chr,
                                                List<CandidateSnp> snps) throws IOException {
        // Build position lookup from candidate SNPs
        Map<String, CandidateSnp> idToSnp = new HashMap<>();
        for (CandidateSnp s : snps) idToSnp.put(s.id, s);

        List<ClumpInterval> clumps = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(clumpedFile.toFile()))) {
            String line;
            boolean headerPassed = false;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (!headerPassed) {
                    if (line.startsWith("CHR") || line.startsWith("chr")) headerPassed = true;
                    continue;
                }
                // PLINK .clumped format: CHR  F  SNP  BP  P  TOTAL  NSIG  S05  S01  S001  S0001  SP2
                String[] f = line.split("\\s+");
                if (f.length < 5) continue;

                String leadSnpId = f[2];
                long leadPos;
                double leadP;
                try {
                    leadPos = Long.parseLong(f[3]);
                    leadP = Double.parseDouble(f[4]);
                } catch (NumberFormatException e) { continue; }

                // Parse SP2 column (member SNPs) to get interval extent
                long minPos = leadPos, maxPos = leadPos;
                int memberCount = 1;
                if (f.length >= 12) {
                    String sp2 = f[11];
                    if (!sp2.equals("NONE")) {
                        String[] members = sp2.replace("(1)", "").replace("(2)", "").split(",");
                        for (String m : members) {
                            m = m.trim();
                            if (m.isEmpty()) continue;
                            CandidateSnp ms = idToSnp.get(m);
                            if (ms != null) {
                                minPos = Math.min(minPos, ms.pos);
                                maxPos = Math.max(maxPos, ms.pos);
                                memberCount++;
                            }
                        }
                    }
                }

                ClumpInterval ci = new ClumpInterval();
                ci.chr = chr;
                ci.start = minPos;
                ci.end = maxPos;
                ci.leadSnp = leadSnpId;
                ci.leadP = leadP;
                ci.memberCount = memberCount;
                clumps.add(ci);
            }
        }
        return clumps;
    }

    // ── Merge clump intervals ──────────────────────────────────────────

    static List<IdentifiedLocus> mergeClumpIntervals(List<ClumpInterval> clumps,
                                                     int mergeDistanceBp, int minSnpsPerLocus) {
        // Group by chromosome
        Map<String, List<ClumpInterval>> byChr = new TreeMap<>(LociIdentifier::chrCompare);
        for (ClumpInterval ci : clumps)
            byChr.computeIfAbsent(ci.chr, k -> new ArrayList<>()).add(ci);

        List<IdentifiedLocus> loci = new ArrayList<>();
        for (Map.Entry<String, List<ClumpInterval>> entry : byChr.entrySet()) {
            String chr = entry.getKey();
            List<ClumpInterval> chrClumps = entry.getValue();
            chrClumps.sort(Comparator.comparingLong(c -> c.start));

            long curStart = chrClumps.get(0).start;
            long curEnd = chrClumps.get(0).end;
            String bestSnp = chrClumps.get(0).leadSnp;
            double bestP = chrClumps.get(0).leadP;
            int totalMembers = chrClumps.get(0).memberCount;

            for (int i = 1; i < chrClumps.size(); i++) {
                ClumpInterval ci = chrClumps.get(i);
                if (ci.start - curEnd <= mergeDistanceBp) {
                    curEnd = Math.max(curEnd, ci.end);
                    totalMembers += ci.memberCount;
                    if (ci.leadP < bestP) { bestP = ci.leadP; bestSnp = ci.leadSnp; }
                } else {
                    if (totalMembers >= minSnpsPerLocus) {
                        IdentifiedLocus l = new IdentifiedLocus();
                        l.chr = chr; l.start = curStart; l.end = curEnd;
                        l.topSnpId = bestSnp; l.topP = bestP; l.nSnps = totalMembers;
                        loci.add(l);
                    }
                    curStart = ci.start; curEnd = ci.end;
                    bestSnp = ci.leadSnp; bestP = ci.leadP; totalMembers = ci.memberCount;
                }
            }
            if (totalMembers >= minSnpsPerLocus) {
                IdentifiedLocus l = new IdentifiedLocus();
                l.chr = chr; l.start = curStart; l.end = curEnd;
                l.topSnpId = bestSnp; l.topP = bestP; l.nSnps = totalMembers;
                loci.add(l);
            }
        }
        return loci;
    }

    // ── Helpers ─────────────────────────────────────────────────────────

    static String findPlink() {
        for (String name : new String[]{"plink", "plink2", "plink.exe", "plink2.exe"}) {
            try {
                Process p = new ProcessBuilder(name, "--version").redirectErrorStream(true).start();
                try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                    while (br.readLine() != null) {}
                }
                if (p.waitFor() == 0) return name;
            } catch (Exception ignored) {}
        }
        // Try common Windows paths
        String[] paths = {
            "C:/Users/alsammana/Documents/DataResources/plink/plink.exe",
            "plink/plink.exe", "tools/plink.exe"
        };
        for (String path : paths) {
            if (new File(path).exists()) return path;
        }
        return null;
    }

    public static void writeLociFiles(List<IdentifiedLocus> loci, String lociPath,
                                      String detailPath) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(lociPath)))) {
            pw.println("meta_chr\tmeta_start\tmeta_end");
            for (IdentifiedLocus l : loci) pw.printf("%s\t%d\t%d%n", l.chr, l.start, l.end);
        }
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(detailPath)))) {
            pw.println(IdentifiedLocus.tsvHeader());
            for (IdentifiedLocus l : loci) pw.println(l.toTsv());
        }
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }

    static int chrToInt(String chr) {
        try { return Integer.parseInt(chr); }
        catch (NumberFormatException e) {
            switch (chr.toUpperCase()) {
                case "X": return 23; case "Y": return 24; case "MT": case "M": return 25;
                default: return 99;
            }
        }
    }

    private static int chrCompare(String a, String b) {
        return Integer.compare(chrToInt(a), chrToInt(b));
    }
}
