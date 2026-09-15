
import java.io.*;
import java.util.*;

/**
 * Step 1: Extract GWAS summary rows for a single locus from the full GWAS file,
 * and extract the matching region from the ancestry reference panel via PLINK.
 *
 * Output:
 *   base/locus_gwas.tsv       — All GWAS rows within the locus padded boundaries
 *   base/ref_region.bed/bim/fam — PLINK bfile subset for the same region
 *   base/base.manifest.json   — Content hashes of inputs
 */
public class LocusGwasExtractor {

    public static class Result {
        public boolean ok;
        public String error;
        public int gwasSnpCount;
        public int refVariantCount;
        public String inputHash;
    }

    /**
     * @return current input hash for cache checking (without running the step)
     */
    public static String computeInputHash(Config config, Locus locus) throws IOException {
        List<String> parts = new ArrayList<>();
        parts.add(ContentHasher.hashFile(new File(config.gwasFile)));
        parts.add(locus.id);
        parts.add(locus.chr + ":" + locus.paddedStart + "-" + locus.paddedEnd);
        parts.add(config.refPanelPath);
        parts.add(config.colChr + "|" + config.colPos + "|" + config.colPvalue +
                   "|" + config.colEa + "|" + config.colNea +
                   "|" + config.colBeta + "|" + config.colSe + "|" + config.colN +
                   "|" + config.colRsid + "|" + config.colVarid +
                   "|" + config.colMaf + "|" + config.colInfo + "|" + config.colOr);
        return ContentHasher.combineHashes(parts);
    }

    public static Result run(Config config, Locus locus, File baseDir) throws IOException {
        Result result = new Result();
        baseDir.mkdirs();

        String inputHash = computeInputHash(config, locus);
        result.inputHash = inputHash;

        // Check cache
        StepManifest existing = StepManifest.read(baseDir, "base");
        if (existing != null && existing.isCurrent(inputHash)) {
            File gwasOut = new File(baseDir, "locus_gwas.tsv");
            File bimOut = new File(baseDir, "ref_region.bim");
            if (gwasOut.exists() && (config.refPanelPath.isEmpty() || bimOut.exists())) {
                result.ok = true;
                System.out.printf("[BaseStep] Locus %s: base step cache hit%n", locus.id);
                return result;
            }
        }

        // ── Extract GWAS rows ──
        File gwasOut = new File(baseDir, "locus_gwas.tsv");
        int snpCount = extractGwasRows(config, locus, gwasOut);
        result.gwasSnpCount = snpCount;
        if (snpCount == 0) {
            result.error = "No GWAS SNPs found in locus region";
            return result;
        }

        // ── Extract ref panel subset ──
        int refCount = 0;
        if (!config.refPanelPath.isEmpty()) {
            String plinkBin = PlinkSubsetter.findPlink(config);
            if (plinkBin == null) {
                result.error = "PLINK binary not found";
                return result;
            }
            refCount = extractRefRegion(plinkBin, config.refPanelPath, locus, baseDir);
            result.refVariantCount = refCount;
            if (refCount == 0) {
                result.error = "No ref panel variants in locus region";
                return result;
            }
        }

        // Write manifest
        StepManifest manifest = new StepManifest();
        manifest.stepName = "base";
        manifest.inputHash = inputHash;
        manifest.refPanelPopulation = config.refPanelPopulation;
        manifest.builtAt = System.currentTimeMillis();
        manifest.inputs.put("gwas_file", config.gwasFile);
        manifest.inputs.put("gwas_file_hash", ContentHasher.hashFile(new File(config.gwasFile)));
        manifest.inputs.put("locus_id", locus.id);
        manifest.inputs.put("region", locus.chr + ":" + locus.paddedStart + "-" + locus.paddedEnd);
        manifest.inputs.put("ref_panel", config.refPanelPath);
        manifest.outputs.put("gwas_snps", String.valueOf(snpCount));
        manifest.outputs.put("ref_variants", String.valueOf(refCount));
        manifest.write(baseDir);

        result.ok = true;
        System.out.printf("[BaseStep] Locus %s: extracted %d GWAS SNPs, %d ref variants%n",
            locus.id, snpCount, refCount);
        return result;
    }

    private static int extractGwasRows(Config config, Locus locus, File outFile)
            throws IOException {
        int count = 0, derivedBetaSe = 0, recomputedPval = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(config.gwasFile), 1024 * 1024);
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile)))) {

            String headerLine = br.readLine();
            if (headerLine == null) return 0;
            String[] header = headerLine.trim().split("\t");

            int iChr  = colIdx(header, config.colChr);
            int iPos  = colIdx(header, config.colPos);
            int iPval = colIdx(header, config.colPvalue);
            int iEa   = colIdx(header, config.colEa);
            int iNea  = colIdx(header, config.colNea);
            if (iChr < 0 || iPos < 0 || iPval < 0 || iEa < 0 || iNea < 0) return 0;

            int iRsid  = config.colRsid.isEmpty()  ? -1 : colIdx(header, config.colRsid);
            int iVarid = config.colVarid.isEmpty()  ? -1 : colIdx(header, config.colVarid);
            int iBeta  = config.colBeta.isEmpty()   ? -1 : colIdx(header, config.colBeta);
            int iOr    = config.colOr.isEmpty()     ? -1 : colIdx(header, config.colOr);
            int iSe    = config.colSe.isEmpty()     ? -1 : colIdx(header, config.colSe);
            int iN     = config.colN.isEmpty()      ? -1 : colIdx(header, config.colN);
            int iMaf   = config.colMaf.isEmpty()    ? -1 : colIdx(header, config.colMaf);
            int iInfo  = config.colInfo.isEmpty()    ? -1 : colIdx(header, config.colInfo);

            // Write standardized header
            pw.println("snp_id\tchr\tpos\tea\tnea\tpvalue\tbeta\tse\tor\tn\tmaf\tinfo\trsid\tvarid");

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length <= Math.max(iChr, Math.max(iPos, Math.max(iPval, Math.max(iEa, iNea)))))
                    continue;

                String chr = f[iChr].trim().replaceFirst("^chr", "");
                if (!chr.equals(locus.chr)) continue;

                long pos;
                try { pos = Long.parseLong(f[iPos].trim()); }
                catch (NumberFormatException e) { continue; }

                if (pos < locus.paddedStart || pos > locus.paddedEnd) continue;

                String ea   = f[iEa].trim().toUpperCase();
                String nea  = f[iNea].trim().toUpperCase();
                String pval = f[iPval].trim();

                String rsid  = safeGet(f, iRsid, "");
                String varid = safeGet(f, iVarid, "");
                if (varid.isEmpty()) varid = chr + ":" + pos;

                // Stable SNP ID: chr:pos:sorted(alleles)
                String a1 = ea.compareTo(nea) <= 0 ? ea : nea;
                String a2 = ea.compareTo(nea) <= 0 ? nea : ea;
                String snpId = chr + ":" + pos + ":" + a1 + ":" + a2;

                String beta = safeGet(f, iBeta, "NA");
                String se   = safeGet(f, iSe, "NA");
                String or_  = safeGet(f, iOr, "NA");
                String n    = safeGet(f, iN, "NA");
                String maf  = safeGet(f, iMaf, "NA");
                String info = safeGet(f, iInfo, "NA");

                // Derive missing beta/SE from OR + p-value (common for older/externally-sourced GWAS
                // files that report only OR) — never overwrites a beta/SE the file already provides.
                if ((beta.equals("NA") || se.equals("NA")) && !or_.equals("NA")) {
                    try {
                        double orVal = Double.parseDouble(or_);
                        double pVal = Double.parseDouble(pval);
                        double[] derived = StatsUtil.deriveBetaSeFromOrP(orVal, pVal);
                        if (derived != null) {
                            beta = String.valueOf(derived[0]);
                            se = String.valueOf(derived[1]);
                            derivedBetaSe++;
                        }
                    } catch (NumberFormatException ignored) {}
                }

                // A literal p=0 (common when a source tool's own p-value underflowed at export time)
                // makes -log10(p) infinite; recompute it from beta/se via the survival function
                // directly rather than trusting a floored/zero value, when both are available.
                if (!beta.equals("NA") && !se.equals("NA")) {
                    double pAsWritten;
                    boolean pLooksFloored;
                    try {
                        pAsWritten = Double.parseDouble(pval);
                        pLooksFloored = !(pAsWritten > 0) || !Double.isFinite(pAsWritten);
                    } catch (NumberFormatException e) { pLooksFloored = true; }
                    if (pLooksFloored) {
                        try {
                            double z = Double.parseDouble(beta) / Double.parseDouble(se);
                            if (Double.isFinite(z) && z != 0) {
                                double recomputed = StatsUtil.zToP(z);
                                if (recomputed > 0) { pval = String.valueOf(recomputed); recomputedPval++; }
                            }
                        } catch (NumberFormatException ignored) {}
                    }
                }

                pw.printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                    snpId, chr, pos, ea, nea, pval, beta, se, or_, n, maf, info, rsid, varid);
                count++;
            }
        }
        if (derivedBetaSe > 0 || recomputedPval > 0) {
            System.out.printf("[BaseStep] Locus %s: derived beta/SE from OR+p for %d SNPs, recomputed %d floored p-values%n",
                locus.id, derivedBetaSe, recomputedPval);
        }
        return count;
    }

    private static int extractRefRegion(String plinkBin, String refPanelPath,
                                         Locus locus, File baseDir) throws IOException {
        String prefix = new File(baseDir, "ref_region").getAbsolutePath();

        List<String> cmd = Arrays.asList(
            plinkBin,
            "--bfile",   refPanelPath,
            "--chr",     locus.chr,
            "--from-bp", String.valueOf(locus.paddedStart),
            "--to-bp",   String.valueOf(locus.paddedEnd),
            "--make-bed",
            "--out",     prefix,
            "--silent"
        );

        try {
            Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                while (br.readLine() != null) { /* drain */ }
            }
            proc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return 0;
        }

        // Count variants from .bim
        File bimFile = new File(prefix + ".bim");
        if (!bimFile.exists()) return 0;

        int count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
            while (br.readLine() != null) count++;
        }

        // Cleanup PLINK temp files
        new File(prefix + ".log").delete();
        new File(prefix + ".nosex").delete();
        new File(prefix + ".noperm").delete();

        return count;
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }

    private static String safeGet(String[] f, int idx, String fallback) {
        if (idx < 0 || idx >= f.length) return fallback;
        String v = f[idx].trim();
        return (v.isEmpty() || v.equals(".") || v.equalsIgnoreCase("NA")) ? fallback : v;
    }
}
