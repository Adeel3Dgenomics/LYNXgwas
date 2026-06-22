
import java.io.*;
import java.util.*;

/**
 * Step 3: Harmonize GWAS so effect alleles align to the reference panel orientation.
 * Emits allele_alignment.tsv recording kept/flipped/dropped per SNP.
 *
 * "Harmonized" means: the effect allele in the GWAS matches the alt (A1) allele
 * in the reference panel. If GWAS effect allele matches the ref allele instead,
 * the beta/OR are flipped.
 */
public class AlleleHarmonizer {

    public static class Result {
        public boolean ok;
        public String error;
        public int kept;
        public int flipped;
        public int complemented;
        public int dropped;
        public String inputHash;
    }

    public static String computeInputHash(File matchedDir) throws IOException {
        return ContentHasher.combineHashes(
            ContentHasher.hashFile(new File(matchedDir, "matched_gwas.tsv")),
            ContentHasher.hashFile(new File(matchedDir, "matched_ref.bim"))
        );
    }

    public static Result run(File matchedDir, File harmonizedDir) throws IOException {
        Result result = new Result();
        harmonizedDir.mkdirs();

        String inputHash = computeInputHash(matchedDir);
        result.inputHash = inputHash;

        StepManifest existing = StepManifest.read(harmonizedDir, "harmonized");
        if (existing != null && existing.isCurrent(inputHash)) {
            if (new File(harmonizedDir, "harmonized_gwas.tsv").exists() &&
                new File(harmonizedDir, "allele_alignment.tsv").exists()) {
                result.ok = true;
                return result;
            }
        }

        // ── Read matched ref .bim for allele orientation ──
        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (!bimFile.exists()) {
            result.error = "matched_ref.bim not found";
            return result;
        }

        Map<String, String[]> refAlleles = new LinkedHashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;
                String chr = f[0].replaceFirst("^chr", "");
                String pos = f[3].trim();
                String a1 = f[4].trim().toUpperCase();
                String a2 = f[5].trim().toUpperCase();
                refAlleles.put(chr + ":" + pos, new String[]{a1, a2});
            }
        }

        // ── Read matched GWAS and harmonize ──
        File gwasIn = new File(matchedDir, "matched_gwas.tsv");
        File gwasOut = new File(harmonizedDir, "harmonized_gwas.tsv");
        File alignOut = new File(harmonizedDir, "allele_alignment.tsv");

        int kept = 0, flipped = 0, complemented = 0, dropped = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(gwasIn));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(gwasOut)));
             PrintWriter al = new PrintWriter(new BufferedWriter(new FileWriter(alignOut)))) {

            String header = br.readLine();
            if (header == null) { result.error = "Empty matched GWAS"; return result; }

            // Header columns: snp_id chr pos ea nea pvalue beta se or n maf info rsid varid
            pw.println(header);
            al.println("snp_id\tchr\tpos\tgwas_ea\tgwas_nea\tref_a1\tref_a2\taction\tharmonized_ea\tharmonized_nea\tbeta_sign");

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;

                String snpId = f[0];
                String chr   = f[1];
                String pos   = f[2];
                String ea    = f[3].toUpperCase();
                String nea   = f[4].toUpperCase();
                String pval  = f[5];
                String beta  = f.length > 6  ? f[6]  : "NA";
                String se    = f.length > 7  ? f[7]  : "NA";
                String or_   = f.length > 8  ? f[8]  : "NA";
                String n     = f.length > 9  ? f[9]  : "NA";
                String maf   = f.length > 10 ? f[10] : "NA";
                String info  = f.length > 11 ? f[11] : "NA";
                String rsid  = f.length > 12 ? f[12] : "";
                String varid = f.length > 13 ? f[13] : "";

                String key = chr + ":" + pos;
                String[] ref = refAlleles.get(key);
                if (ref == null) { dropped++; continue; }
                String refA1 = ref[0], refA2 = ref[1];

                SnpMatcher.AlleleMatch match = SnpMatcher.matchAlleles(ea, nea, refA1, refA2);
                if (match == SnpMatcher.AlleleMatch.NONE) {
                    al.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\tdropped\t\t\t%n",
                        snpId, chr, pos, ea, nea, refA1, refA2);
                    dropped++;
                    continue;
                }

                String harmEa = ea, harmNea = nea;
                String action = "kept";
                int betaSign = 1;

                // If GWAS ea matches ref A2 (not A1), flip
                boolean needFlip = false;
                if (match == SnpMatcher.AlleleMatch.FORWARD || match == SnpMatcher.AlleleMatch.COMPLEMENT) {
                    if (!ea.equals(refA1) && !SnpMatcher.complement(ea).equals(refA1)) {
                        needFlip = true;
                    }
                } else {
                    // REVERSE or REVERSE_COMPLEMENT: ea matched refA2
                    needFlip = true;
                }

                if (match == SnpMatcher.AlleleMatch.COMPLEMENT || match == SnpMatcher.AlleleMatch.REVERSE_COMPLEMENT) {
                    harmEa = SnpMatcher.complement(ea);
                    harmNea = SnpMatcher.complement(nea);
                    action = "complemented";
                    complemented++;
                }

                if (needFlip) {
                    // Swap ea/nea and negate beta
                    String tmp = harmEa; harmEa = harmNea; harmNea = tmp;
                    betaSign = -1;
                    action = action.equals("complemented") ? "complement_flipped" : "flipped";
                    flipped++;

                    beta = flipValue(beta);
                    or_ = flipOr(or_);
                } else {
                    kept++;
                }

                // Write harmonized row
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                    snpId, chr, pos, harmEa, harmNea, pval, beta, se, or_, n, maf, info, rsid, varid);

                al.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d%n",
                    snpId, chr, pos, ea, nea, refA1, refA2, action, harmEa, harmNea, betaSign);
            }
        }

        result.kept = kept;
        result.flipped = flipped;
        result.complemented = complemented;
        result.dropped = dropped;
        result.ok = true;

        StepManifest manifest = new StepManifest();
        manifest.stepName = "harmonized";
        manifest.inputHash = inputHash;
        manifest.builtAt = System.currentTimeMillis();
        manifest.outputs.put("kept", String.valueOf(kept));
        manifest.outputs.put("flipped", String.valueOf(flipped));
        manifest.outputs.put("complemented", String.valueOf(complemented));
        manifest.outputs.put("dropped", String.valueOf(dropped));
        manifest.write(harmonizedDir);

        System.out.printf("[Harmonizer] %d kept, %d flipped, %d complemented, %d dropped%n",
            kept, flipped, complemented, dropped);
        return result;
    }

    private static String flipValue(String val) {
        if (val == null || val.equals("NA") || val.isEmpty()) return val;
        try {
            double d = Double.parseDouble(val);
            return Double.isNaN(d) ? val : String.valueOf(-d);
        } catch (NumberFormatException e) { return val; }
    }

    private static String flipOr(String val) {
        if (val == null || val.equals("NA") || val.isEmpty()) return val;
        try {
            double d = Double.parseDouble(val);
            return (Double.isNaN(d) || d == 0) ? val : String.valueOf(1.0 / d);
        } catch (NumberFormatException e) { return val; }
    }
}
