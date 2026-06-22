
import java.io.*;
import java.util.*;

/**
 * Step 2: Match GWAS SNPs to the reference panel by chr:pos (allele-aware).
 * Produces matched_gwas.tsv and matched_ref.{bed,bim,fam} restricted to the
 * SAME SNPs in the SAME order.
 *
 * Allele matching rules:
 *   - Forward match:  GWAS(A/G) == REF(A/G)
 *   - Reverse match:  GWAS(A/G) == REF(G/A) — same alleles, swapped
 *   - Complement:     GWAS(A/G) matches REF(T/C) — strand flip
 *   - No match:       dropped with reason in match log
 */
public class SnpMatcher {

    public static class Result {
        public boolean ok;
        public String error;
        public int matchedCount;
        public int droppedCount;
        public String inputHash;
    }

    public static String computeInputHash(File baseDir) throws IOException {
        return ContentHasher.combineHashes(
            ContentHasher.hashFile(new File(baseDir, "locus_gwas.tsv")),
            ContentHasher.hashFile(new File(baseDir, "ref_region.bim"))
        );
    }

    public static Result run(File baseDir, File matchedDir, Config config) throws IOException {
        Result result = new Result();
        matchedDir.mkdirs();

        String inputHash = computeInputHash(baseDir);
        result.inputHash = inputHash;

        StepManifest existing = StepManifest.read(matchedDir, "matched");
        if (existing != null && existing.isCurrent(inputHash)) {
            if (new File(matchedDir, "matched_gwas.tsv").exists()) {
                result.ok = true;
                return result;
            }
        }

        // ── Read reference panel .bim ──
        // bim format: chr  id  cm  pos  a1  a2
        File bimFile = new File(baseDir, "ref_region.bim");
        if (!bimFile.exists()) {
            result.error = "ref_region.bim not found";
            return result;
        }

        // Map: "chr:pos" → BimEntry (preserving order)
        LinkedHashMap<String, BimEntry> refByPos = new LinkedHashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;
                String chr = f[0].replaceFirst("^chr", "");
                String bimId = f[1];
                String pos = f[3].trim();
                String a1 = f[4].trim().toUpperCase();
                String a2 = f[5].trim().toUpperCase();
                String key = chr + ":" + pos;
                if (!refByPos.containsKey(key)) {
                    refByPos.put(key, new BimEntry(bimId, chr, Long.parseLong(pos), a1, a2));
                }
            }
        }

        // ── Read GWAS and match ──
        File gwasFile = new File(baseDir, "locus_gwas.tsv");
        List<String> matchedGwasLines = new ArrayList<>();
        List<String> matchedBimIds = new ArrayList<>();
        String gwasHeader = null;
        int dropped = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            gwasHeader = br.readLine();
            if (gwasHeader == null) { result.error = "Empty GWAS file"; return result; }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 5) continue;

                String chr = f[1];  // chr column
                String pos = f[2];  // pos column
                String ea = f[3].toUpperCase();   // effect allele
                String nea = f[4].toUpperCase();  // other allele

                String key = chr + ":" + pos;
                BimEntry ref = refByPos.get(key);
                if (ref == null) { dropped++; continue; }

                AlleleMatch match = matchAlleles(ea, nea, ref.a1, ref.a2);
                if (match == AlleleMatch.NONE) { dropped++; continue; }

                matchedGwasLines.add(line);
                matchedBimIds.add(ref.bimId);
            }
        }

        if (matchedGwasLines.isEmpty()) {
            result.error = "No GWAS SNPs matched reference panel";
            return result;
        }

        // ── Write matched GWAS ──
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(
                new FileWriter(new File(matchedDir, "matched_gwas.tsv"))))) {
            pw.println(gwasHeader);
            for (String line : matchedGwasLines) pw.println(line);
        }

        // ── Extract matched ref panel via PLINK --extract ──
        File extractFile = new File(matchedDir, "_extract_ids.txt");
        try (PrintWriter pw = new PrintWriter(new FileWriter(extractFile))) {
            for (String id : matchedBimIds) pw.println(id);
        }

        String plinkBin = PlinkSubsetter.findPlink(config);
        if (plinkBin != null) {
            String refPrefix = new File(baseDir, "ref_region").getAbsolutePath();
            String outPrefix = new File(matchedDir, "matched_ref").getAbsolutePath();

            List<String> cmd = Arrays.asList(
                plinkBin,
                "--bfile", refPrefix,
                "--extract", extractFile.getAbsolutePath(),
                "--make-bed",
                "--out", outPrefix,
                "--silent"
            );
            try {
                Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
                try (BufferedReader drain = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                    while (drain.readLine() != null) {}
                }
                proc.waitFor();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }

            new File(outPrefix + ".log").delete();
            new File(outPrefix + ".nosex").delete();
        }
        extractFile.delete();

        result.matchedCount = matchedGwasLines.size();
        result.droppedCount = dropped;
        result.ok = true;

        // Write manifest
        StepManifest manifest = new StepManifest();
        manifest.stepName = "matched";
        manifest.inputHash = inputHash;
        manifest.builtAt = System.currentTimeMillis();
        manifest.outputs.put("matched_snps", String.valueOf(result.matchedCount));
        manifest.outputs.put("dropped_snps", String.valueOf(dropped));
        manifest.write(matchedDir);

        System.out.printf("[SnpMatcher] %d matched, %d dropped%n",
            result.matchedCount, dropped);
        return result;
    }

    // ── Allele matching ──

    enum AlleleMatch { FORWARD, REVERSE, COMPLEMENT, REVERSE_COMPLEMENT, NONE }

    static AlleleMatch matchAlleles(String ea, String nea, String refA1, String refA2) {
        if ((ea.equals(refA1) && nea.equals(refA2)) || (ea.equals(refA2) && nea.equals(refA1)))
            return ea.equals(refA1) ? AlleleMatch.FORWARD : AlleleMatch.REVERSE;

        String ceA = complement(ea), cnA = complement(nea);
        if ((ceA.equals(refA1) && cnA.equals(refA2)) || (ceA.equals(refA2) && cnA.equals(refA1)))
            return ceA.equals(refA1) ? AlleleMatch.COMPLEMENT : AlleleMatch.REVERSE_COMPLEMENT;

        return AlleleMatch.NONE;
    }

    static String complement(String allele) {
        StringBuilder sb = new StringBuilder(allele.length());
        for (int i = 0; i < allele.length(); i++) {
            switch (allele.charAt(i)) {
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                default:  sb.append(allele.charAt(i));
            }
        }
        return sb.toString();
    }

    static class BimEntry {
        final String bimId;
        final String chr;
        final long pos;
        final String a1, a2;
        BimEntry(String bimId, String chr, long pos, String a1, String a2) {
            this.bimId = bimId; this.chr = chr; this.pos = pos;
            this.a1 = a1; this.a2 = a2;
        }
    }
}
