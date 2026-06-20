package rsid;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

// TODO: BUILD MISMATCH — verify/liftover build before production use.

public class RsidRecovery {

    private final String snpdbDir;
    private final String build;

    public RsidRecovery(String snpdbDir, String build) {
        this.snpdbDir = snpdbDir;
        this.build    = build;
    }

    /**
     * Build a position index for a genomic region from the dbSNP VCF.
     */
    public Map<Integer, List<DbSnpRecord>> buildPositionIndex(String chr, int start, int end) {
        Map<Integer, List<DbSnpRecord>> posIndex = new LinkedHashMap<>();
        Path vcfPath = Paths.get(snpdbDir, chr + ".vcf.gz");

        if (!Files.exists(vcfPath)) {
            System.err.printf("[RsidRecovery] VCF not found: %s%n", vcfPath);
            return posIndex;
        }

        try (VCFFileReader reader = new VCFFileReader(vcfPath, true)) {
            try (htsjdk.samtools.util.CloseableIterator<VariantContext> it =
                     reader.query(chr, start, end)) {
                while (it.hasNext()) {
                    VariantContext vc = it.next();
                    String rsid = vc.getID();
                    if (rsid == null || rsid.equals(".")) continue;

                    String ref = vc.getReference().getBaseString().toUpperCase();
                    List<String> alts = vc.getAlternateAlleles().stream()
                        .map(a -> a.getBaseString().toUpperCase())
                        .collect(Collectors.toList());

                    Integer dbSnpBuild = null;
                    Object buildAttr = vc.getAttribute("dbSNPBuildID");
                    if (buildAttr != null) {
                        try { dbSnpBuild = Integer.parseInt(buildAttr.toString()); }
                        catch (NumberFormatException ignored) {}
                    }

                    DbSnpRecord rec = new DbSnpRecord(rsid, ref, alts, dbSnpBuild);
                    posIndex.computeIfAbsent(vc.getStart(), k -> new ArrayList<>()).add(rec);
                }
            }
        } catch (Exception e) {
            System.err.printf("[RsidRecovery] Error querying %s:%d-%d: %s%n",
                chr, start, end, e.getMessage());
        }
        return posIndex;
    }

    /**
     * Match a list of GWAS SNPs against a prebuilt position index.
     */
    public List<MatchResult> matchSnps(String chr, List<int[]> snps,
                                       Map<Integer, List<DbSnpRecord>> posIndex) {
        List<MatchResult> results = new ArrayList<>();
        for (int[] snp : snps) {
            // snp = {pos, eaIdx, neaIdx} — but we pass strings directly
            // This overload is for the main entry point below
        }
        return results;
    }

    /**
     * Full recovery: read GWAS file, build index per locus, match, and report.
     */
    public List<MatchResult> recover(String gwasFile, String lociFile,
                                     String colChr, String colPos,
                                     String colEa, String colNea) throws IOException {
        // Parse loci
        List<int[]> loci = new ArrayList<>(); // {chr_int, start, end}
        List<String> lociChrs = new ArrayList<>();
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

        // Parse GWAS file — collect SNPs per locus
        List<String[]> gwasSnps = new ArrayList<>(); // {chr, pos, ea, nea}
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            String[] hcols = header.trim().split("\t");
            int iChr2 = colIdx(hcols, colChr), iPos = colIdx(hcols, colPos);
            int iEa = colIdx(hcols, colEa), iNea = colIdx(hcols, colNea);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t");
                if (f.length <= Math.max(iChr2, Math.max(iPos, Math.max(iEa, iNea)))) continue;
                String chr = f[iChr2].trim();
                String posStr = f[iPos].trim();
                String ea = f[iEa].trim();
                String nea = f[iNea].trim();
                int pos;
                try { pos = Integer.parseInt(posStr); }
                catch (NumberFormatException e) { continue; }

                // Check if this SNP falls in any locus
                for (int li = 0; li < loci.size(); li++) {
                    if (chr.equals(lociChrs.get(li))
                            && pos >= loci.get(li)[0] && pos <= loci.get(li)[1]) {
                        gwasSnps.add(new String[]{chr, posStr, ea, nea});
                        break;
                    }
                }
            }
        }

        System.out.printf("[RsidRecovery] %d GWAS SNPs in %d loci%n", gwasSnps.size(), loci.size());

        // Build position index per chromosome-region and match
        List<MatchResult> results = new ArrayList<>();
        for (int li = 0; li < loci.size(); li++) {
            String chr = lociChrs.get(li);
            int start = loci.get(li)[0], end = loci.get(li)[1];
            System.out.printf("[RsidRecovery] Locus %d: %s:%d-%d ... ", li + 1, chr, start, end);

            Map<Integer, List<DbSnpRecord>> posIndex = buildPositionIndex(chr, start, end);
            System.out.printf("%d dbSNP positions loaded%n", posIndex.size());

            // Match SNPs in this locus
            for (String[] snp : gwasSnps) {
                if (!snp[0].equals(chr)) continue;
                int pos = Integer.parseInt(snp[1]);
                if (pos < start || pos > end) continue;
                results.add(RsidMatcher.match(chr, pos, snp[2], snp[3], posIndex));
            }
        }

        return results;
    }

    /**
     * Match GWAS SNPs directly (for test harness with pre-built data).
     */
    public static List<MatchResult> matchDirect(String chr, String[][] gwasSnps,
                                                Map<Integer, List<DbSnpRecord>> posIndex) {
        List<MatchResult> results = new ArrayList<>();
        for (String[] snp : gwasSnps) {
            int pos = Integer.parseInt(snp[1]);
            results.add(RsidMatcher.match(chr, pos, snp[2], snp[3], posIndex));
        }
        return results;
    }

    /**
     * Print summary statistics.
     */
    public static void printSummary(List<MatchResult> results) {
        int total = results.size();
        int forward = 0, reverse = 0, posMismatch = 0, noPos = 0, multiAllelic = 0, multiVersion = 0;
        for (MatchResult r : results) {
            switch (r.matchReason) {
                case "matched_forward":          forward++;    break;
                case "matched_reverse":          reverse++;    break;
                case "pos_only_allele_mismatch": posMismatch++; break;
                case "no_pos":                   noPos++;      break;
            }
            if (r.nCandidatesAtPos > 1) multiVersion++;
        }
        int matched = forward + reverse;
        System.out.println();
        System.out.printf("Total GWAS SNPs processed:           %d%n", total);
        System.out.printf("Matched (forward + reverse):         %d  (%.1f%%)%n", matched, pct(matched, total));
        System.out.printf("  - forward orientation:             %d%n", forward);
        System.out.printf("  - reverse orientation:             %d%n", reverse);
        System.out.printf("Position found but allele mismatch:  %d  (%.1f%%)%n", posMismatch, pct(posMismatch, total));
        System.out.printf("Position not found in dbSNP:         %d  (%.1f%%)%n", noPos, pct(noPos, total));
        System.out.printf("Positions with >1 rsID version:      %d%n", multiVersion);
        System.out.printf("%nRecovery rate: %.1f%%%n", pct(matched, total));
    }

    private static double pct(int n, int total) {
        return total > 0 ? n * 100.0 / total : 0;
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }

    // NCBI API fallback — STUB ONLY (do not implement)
    // TODO (later): for SNPs still unmatched after local dbSNP matching,
    // query NCBI dbSNP via API to recover rsID.
    // NOT IMPLEMENTED YET — API is rate-limited and error-prone; add only
    // after local matching is verified. In the real pipeline this must be
    // restricted to LEAD SNPs only, to stay within rate limits.
    public static List<MatchResult> ncbiApiFallback(List<MatchResult> unmatched, String build) {
        return unmatched; // stub: returns input unchanged
    }
}
