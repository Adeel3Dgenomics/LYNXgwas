
import java.io.*;
import java.nio.file.Files;
import java.util.*;

/**
 * MAGMA (de Leeuw et al. 2015, https://ctg.cncr.nl/software/magma) gene-based association and
 * gene-set analysis adapter.
 *
 * Three MAGMA invocations, matching MAGMA's own three-step design:
 *   1. annotate()       — SNP-to-gene mapping (--annotate), producing a .genes.annot file.
 *   2. geneAnalysis()   — the gene-based test itself (--gene-annot + --pval + --bfile),
 *                         producing .genes.out (per-gene P) and .genes.raw (for step 3).
 *   3. geneSetAnalysis() — competitive gene-set enrichment (--gene-results + --set-annot),
 *                         producing .gsa.out.
 *
 * Inputs are built from data this app already has for a locus, following the same conventions
 * as CojoAdapter/FinemapAdapter:
 *   - matchedDir/matched_ref.bim              -> MAGMA's SNP location file (SNP, CHR, BP) and
 *                                                the --bfile reference panel for LD.
 *   - harmonizedDir/harmonized_gwas.tsv       -> MAGMA's .pval file (SNP, P[, N]).
 *   - the app's existing Gene/GffParser model -> MAGMA's gene location file (GENE, CHR, START,
 *                                                STOP, STRAND), reused as-is rather than
 *                                                inventing a new gene record type.
 *
 * Unlike CojoAdapter/FinemapAdapter, MAGMA's gene-based test consumes raw P-values (a one-sided
 * SNP-wise mean chi-square/Z statistic), not signed beta+frequency — so there is no allele
 * orientation to verify here (a P-value has no direction to flip). The one correctness check
 * that *does* carry over from CojoAdapter's Fix 2 ("never fall back to the GWAS's own SNP id")
 * is requiring a reference-panel ID match by chr:pos before writing a row to the .pval file,
 * since MAGMA's --pval SNP column must agree with the IDs used in --snp-loc (the ref panel's
 * own BIM IDs); SNPs with no ref panel match are dropped and counted, exactly like CojoAdapter
 * drops SNPs with no ref ID.
 */
public class MagmaAdapter {

    /** Default annotation window (kb upstream+downstream of each gene). Matches this app's
     *  Config.locusPadding convention (200_000 bp = 200 kb padding around a locus), so the
     *  "how far from a gene is still potentially relevant" assumption is consistent with the
     *  rest of the pipeline. Callers may pass any other value to annotate(). */
    public static final long DEFAULT_WINDOW_KB = 200;

    // ── Result records ──────────────────────────────────────────────────────────────────

    /** One row of MAGMA's .genes.out (gene-based test result). */
    public static class GeneResult {
        public String gene;
        public String chr;
        public long start;
        public long stop;
        public int nsnps;
        public int nparam;
        public int n;
        public double zstat;
        public double p;
    }

    /** One row of MAGMA's .gsa.out (competitive gene-set test result). */
    public static class GeneSetResult {
        public String setId;
        public int nGenes;
        public double beta;
        public double betaStd;
        public double se;
        public double p;
    }

    public static class AnnotateResult {
        public boolean ok;
        public String error;
        public File snpLocFile;
        public File geneLocFile;
        public File annotFile;   // <prefix>.genes.annot — required input to geneAnalysis()
        public int nSnps;
        public int nGenes;
        public String logTail;
    }

    public static class GeneAnalysisResult {
        public boolean ok;
        public String error;
        public List<GeneResult> genes = new ArrayList<>();
        public File genesOutFile;   // .genes.out
        public File genesRawFile;   // .genes.raw — required input to geneSetAnalysis()
        public String manifestPath;
        public int snpsWritten;
        public int snpsDroppedNoRefId;
        public int snpsDroppedStats;
        public String logTail;
    }

    public static class GeneSetAnalysisResult {
        public boolean ok;
        public String error;
        public List<GeneSetResult> sets = new ArrayList<>();
        public File gsaOutFile;
        public String manifestPath;
        public String logTail;
    }

    // ── Binary resolution ───────────────────────────────────────────────────────────────

    /**
     * Resolves the MAGMA binary: appRoot/bin/magma(.exe) first, then PATH. Follows the same
     * config-value -> remembered-path -> PATH order used elsewhere in this app for external
     * tools (see PlinkSubsetter.findPlink() for the bin/ check; PATH search added here since
     * MAGMA is not normally bundled with a PLINK-style local copy).
     */
    public static String findMagmaBinary(File appRoot) throws IOException {
        boolean windows = System.getProperty("os.name", "").toLowerCase().contains("win");
        String exeName = windows ? "magma.exe" : "magma";

        File inBin = new File(appRoot, "bin" + File.separator + exeName);
        if (inBin.exists()) return inBin.getAbsolutePath();

        String path = System.getenv("PATH");
        if (path != null) {
            for (String dir : path.split(File.pathSeparator)) {
                if (dir.isEmpty()) continue;
                File candidate = new File(dir, exeName);
                if (candidate.exists()) return candidate.getAbsolutePath();
            }
        }

        throw new IOException("MAGMA binary not found at: " + inBin.getAbsolutePath()
            + " (or on PATH). Place magma.exe in the bin/ folder.");
    }

    /**
     * Verifies an already-resolved MAGMA binary path exists before launching a process with it,
     * exactly like CojoAdapter verifies gctaBin up front (a clear, specific error beats a native
     * ProcessBuilder "CreateProcess error=2" message). Does not itself search PATH — that is
     * findMagmaBinary()'s job; this only re-checks whatever path the caller already resolved.
     */
    private static void verifyMagmaBinary(String magmaBin) throws IOException {
        File f = new File(magmaBin);
        if (!f.isAbsolute()) f = f.getAbsoluteFile();
        if (!f.exists())
            throw new IOException("MAGMA binary not found at: " + f.getAbsolutePath()
                + ". Place magma.exe in the bin/ folder, or resolve one via findMagmaBinary().");
    }

    // ── Step 1: --annotate ──────────────────────────────────────────────────────────────

    /**
     * Builds MAGMA's SNP location and gene location input files for the given locus and runs
     * --annotate to produce the .genes.annot SNP-to-gene mapping.
     *
     * @param matchedDir matched_ref.{bed,bim,fam} directory (already subset+matched to this locus)
     * @param gff        already-loaded gene annotation (GffParser.parse(config))
     * @param runDir     output directory for this MAGMA run
     * @param chr        locus chromosome (with or without "chr" prefix)
     * @param locusStart locus start position (bp)
     * @param locusEnd   locus end position (bp)
     * @param windowKb   up/downstream window MAGMA adds around each gene (see DEFAULT_WINDOW_KB)
     * @param magmaBin   resolved path to the magma binary (see findMagmaBinary)
     */
    public static AnnotateResult annotate(File matchedDir, GffParser gff, File runDir,
                                           String chr, long locusStart, long locusEnd,
                                           long windowKb, String magmaBin) throws IOException {
        runDir.mkdirs();
        verifyMagmaBinary(magmaBin);
        String chrNorm = chr.replaceFirst("^chr", "");

        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (!bimFile.exists())
            throw new IOException("matched_ref.bim not found in " + matchedDir.getAbsolutePath()
                + ". Run the base pipeline first to build the matched reference panel.");

        // ── SNP location file: SNP CHR BP (no header — MAGMA's --snp-loc format) ──
        File snpLoc = new File(runDir, "magma_snploc.txt");
        int nSnps = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(bimFile));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(snpLoc)))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 4) continue;
                String snpChr = f[0].replaceFirst("^chr", "");
                pw.printf("%s\t%s\t%s%n", f[1], snpChr, f[3].trim());
                nSnps++;
            }
        }
        if (nSnps == 0) throw new IOException("No SNPs found in " + bimFile.getAbsolutePath());

        // ── Gene location file: GENE CHR START STOP STRAND (no header) ──
        // Unpadded gene boundaries — MAGMA's own --annotate window=<kb> flag applies the
        // up/downstream padding, so we hand it the genes' true coordinates here.
        List<Gene> genes = gff.overlapping(chrNorm, locusStart, locusEnd);
        File geneLoc = new File(runDir, "magma_geneloc.txt");
        int nGenes = 0;
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(geneLoc)))) {
            for (Gene g : genes) {
                String strand = (g.strand != null && !g.strand.isEmpty()) ? g.strand : "+";
                pw.printf("%s\t%s\t%d\t%d\t%s%n", g.geneId, g.chr, g.start, g.end, strand);
                nGenes++;
            }
        }

        AnnotateResult result = new AnnotateResult();
        result.snpLocFile = snpLoc;
        result.geneLocFile = geneLoc;
        result.nSnps = nSnps;
        result.nGenes = nGenes;

        if (nGenes == 0) {
            result.ok = false;
            result.error = "No genes overlapping " + chrNorm + ":" + locusStart + "-" + locusEnd
                + " in the loaded gene annotation. Nothing to annotate.";
            return result;
        }

        File prefix = new File(runDir, "magma_step1");
        List<String> cmd = Arrays.asList(magmaBin,
            "--annotate", "window=" + windowKb,
            "--snp-loc", snpLoc.getAbsolutePath(),
            "--gene-loc", geneLoc.getAbsolutePath(),
            "--out", prefix.getAbsolutePath());

        File logFile = new File(runDir, "magma_annotate.log");
        int exit = runProcess(cmd, runDir, logFile);
        result.logTail = tailOf(logFile, 60);

        File annotFile = new File(prefix.getAbsolutePath() + ".genes.annot");
        if (exit != 0 || !annotFile.exists()) {
            result.ok = false;
            result.error = "MAGMA --annotate exited with code " + exit + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }
        result.annotFile = annotFile;
        result.ok = true;
        return result;
    }

    // ── Step 2: gene-based test ─────────────────────────────────────────────────────────

    /**
     * Runs MAGMA's gene-based association test: --bfile <ref> --pval <file> [N=<n> | use=...
     * ncol=N] --gene-annot <annot from annotate()> --out <prefix>, then parses .genes.out.
     */
    public static GeneAnalysisResult geneAnalysis(File harmonizedDir, File matchedDir,
                                                   AnnotateResult annotateResult, File runDir,
                                                   int sampleN, String magmaBin) throws IOException {
        if (annotateResult == null || !annotateResult.ok || annotateResult.annotFile == null)
            throw new IOException("Gene annotation (annotate()) must complete successfully before geneAnalysis().");
        runDir.mkdirs();
        verifyMagmaBinary(magmaBin);

        PvalBuildResult pvalBuild = buildPvalFile(harmonizedDir, matchedDir, runDir, sampleN);

        String bfilePrefix = new File(matchedDir, "matched_ref").getAbsolutePath();
        File outPrefix = new File(runDir, "magma_step2");
        List<String> cmd = new ArrayList<>(Arrays.asList(magmaBin,
            "--bfile", bfilePrefix,
            "--pval", pvalBuild.pvalFile.getAbsolutePath()));
        if (pvalBuild.hasPerSnpN) {
            cmd.add("use=SNP,P");
            cmd.add("ncol=N");
        } else {
            cmd.add("N=" + sampleN);
        }
        cmd.add("--gene-annot");
        cmd.add(annotateResult.annotFile.getAbsolutePath());
        cmd.add("--out");
        cmd.add(outPrefix.getAbsolutePath());

        File logFile = new File(runDir, "magma_gene_analysis.log");
        GeneAnalysisResult result = new GeneAnalysisResult();
        result.snpsWritten = pvalBuild.written;
        result.snpsDroppedNoRefId = pvalBuild.droppedNoRefId;
        result.snpsDroppedStats = pvalBuild.droppedStats;

        int exit = runProcess(cmd, runDir, logFile);
        result.logTail = tailOf(logFile, 60);

        File genesOut = new File(outPrefix.getAbsolutePath() + ".genes.out");
        File genesRaw = new File(outPrefix.getAbsolutePath() + ".genes.raw");
        if (exit != 0 || !genesOut.exists()) {
            result.ok = false;
            result.error = "MAGMA gene analysis exited with code " + exit + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }
        result.genesOutFile = genesOut;
        if (genesRaw.exists()) result.genesRawFile = genesRaw;
        result.genes = parseGenesOut(genesOut);
        result.ok = true;

        File manifestFile = new File(runDir, "magma_genes.manifest.json");
        writeGeneManifest(manifestFile, result.genes.size(), pvalBuild.written, pvalBuild.droppedNoRefId,
            pvalBuild.droppedStats, pvalBuild.hasPerSnpN, sampleN);
        result.manifestPath = manifestFile.getAbsolutePath();

        return result;
    }

    /** Holds the outcome of writing MAGMA's .pval input file — split out of geneAnalysis() so the
     *  ref-panel-ID-matching/N-handling logic can be exercised by tests without launching MAGMA. */
    static class PvalBuildResult {
        File pvalFile;
        int written;
        int droppedNoRefId;
        int droppedStats;
        boolean hasPerSnpN;
    }

    /**
     * Writes MAGMA's .pval file (SNP P [N]) from harmonized_gwas.tsv, matching each row to the
     * reference panel's own BIM SNP ID by chr:pos (never falling back to the GWAS's own SNP id —
     * same correctness reason as CojoAdapter's Fix 2, since MAGMA's --pval SNP column must agree
     * with the IDs used in --snp-loc). Rows with no ref panel match, or a missing/invalid P
     * (and, when using per-SNP N, no usable N), are dropped and counted rather than written.
     */
    static PvalBuildResult buildPvalFile(File harmonizedDir, File matchedDir, File runDir, int sampleN)
            throws IOException {
        runDir.mkdirs();
        File gwasFile = new File(harmonizedDir, "harmonized_gwas.tsv");
        if (!gwasFile.exists())
            throw new IOException("harmonized_gwas.tsv not found in " + harmonizedDir.getAbsolutePath());

        // Ref panel chr:pos -> BIM SNP id.
        Map<String, String> posToRefId = new LinkedHashMap<>();
        File bimFile = new File(matchedDir, "matched_ref.bim");
        if (bimFile.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(bimFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] f = line.split("\t", -1);
                    if (f.length < 4) continue;
                    String bimChr = f[0].replaceFirst("^chr", "");
                    posToRefId.put(bimChr + ":" + f[3].trim(), f[1]);
                }
            }
        }

        // First pass: does the GWAS carry a usable per-SNP N column?
        boolean hasPerSnpN = false;
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            String header = br.readLine();
            if (header == null) throw new IOException("harmonized_gwas.tsv is empty.");
            int iN = findColumn(header, "n");
            if (iN >= 0) {
                String line = br.readLine();
                if (line != null) {
                    String[] f = line.split("\t", -1);
                    if (iN < f.length) {
                        String nVal = f[iN].trim();
                        if (!nVal.isEmpty() && !nVal.equals("NA")) {
                            try { if ((int) Double.parseDouble(nVal) > 0) hasPerSnpN = true; }
                            catch (NumberFormatException ignored) {}
                        }
                    }
                }
            }
        }

        if (sampleN <= 0 && !hasPerSnpN)
            throw new IOException("Sample size (N) is required for MAGMA gene analysis. "
                + "Set it in the project configuration or ensure the GWAS has a per-SNP N column.");

        PvalBuildResult result = new PvalBuildResult();
        result.hasPerSnpN = hasPerSnpN;
        result.pvalFile = new File(runDir, "magma_input.pval");

        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile));
             PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(result.pvalFile)))) {
            String header = br.readLine();
            int iN = findColumn(header, "n");

            pw.println(hasPerSnpN ? "SNP\tP\tN" : "SNP\tP");

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length < 6) continue;
                String chr = f[1], pos = f[2], pval = f[5];
                if (pval.equals("NA")) { result.droppedStats++; continue; }
                double p;
                try { p = Double.parseDouble(pval); } catch (NumberFormatException e) { result.droppedStats++; continue; }
                if (Double.isNaN(p) || p <= 0 || p > 1) { result.droppedStats++; continue; }

                String refId = posToRefId.get(chr.replaceFirst("^chr", "") + ":" + pos);
                if (refId == null) { result.droppedNoRefId++; continue; }

                if (hasPerSnpN) {
                    int snpN = sampleN > 0 ? sampleN : 0;
                    if (iN >= 0 && iN < f.length && !f[iN].trim().isEmpty() && !f[iN].trim().equals("NA")) {
                        try {
                            int parsed = (int) Double.parseDouble(f[iN].trim());
                            if (parsed > 0) snpN = parsed;
                        } catch (NumberFormatException ignored) {}
                    }
                    if (snpN <= 0) { result.droppedStats++; continue; }
                    pw.printf("%s\t%.6g\t%d%n", refId, p, snpN);
                } else {
                    pw.printf("%s\t%.6g%n", refId, p);
                }
                result.written++;
            }
        }

        if (result.written == 0)
            throw new IOException("No usable SNPs (valid P-value + reference panel match) written to .pval file.");

        return result;
    }

    // ── Step 3: gene-set analysis (optional) ───────────────────────────────────────────

    /**
     * Runs MAGMA's competitive gene-set test: --gene-results <genes.raw from geneAnalysis()>
     * --set-annot <setAnnotFile> --out <prefix>, then parses .gsa.out.
     *
     * @param setAnnotFile a MAGMA gene-set annotation file (one set per line: SET_ID GENE GENE ...)
     */
    public static GeneSetAnalysisResult geneSetAnalysis(File genesRawFile, File setAnnotFile,
                                                          File runDir, String magmaBin) throws IOException {
        if (genesRawFile == null || !genesRawFile.exists())
            throw new IOException("Gene-based .genes.raw file not found"
                + (genesRawFile != null ? " at " + genesRawFile.getAbsolutePath() : "")
                + " — run geneAnalysis() first; it is required input for --gene-results.");
        if (setAnnotFile == null || !setAnnotFile.exists())
            throw new IOException("Gene-set annotation file not found"
                + (setAnnotFile != null ? ": " + setAnnotFile.getAbsolutePath() : "."));
        runDir.mkdirs();
        verifyMagmaBinary(magmaBin);

        File outPrefix = new File(runDir, "magma_step3");
        List<String> cmd = Arrays.asList(magmaBin,
            "--gene-results", genesRawFile.getAbsolutePath(),
            "--set-annot", setAnnotFile.getAbsolutePath(),
            "--out", outPrefix.getAbsolutePath());

        File logFile = new File(runDir, "magma_geneset_analysis.log");
        GeneSetAnalysisResult result = new GeneSetAnalysisResult();
        int exit = runProcess(cmd, runDir, logFile);
        result.logTail = tailOf(logFile, 60);

        File gsaOut = new File(outPrefix.getAbsolutePath() + ".gsa.out");
        if (exit != 0 || !gsaOut.exists()) {
            result.ok = false;
            result.error = "MAGMA gene-set analysis exited with code " + exit + ". Check " + logFile.getAbsolutePath()
                + (result.logTail != null ? "\n---\n" + result.logTail : "");
            return result;
        }
        result.gsaOutFile = gsaOut;
        result.sets = parseGsaOut(gsaOut);
        result.ok = true;

        File manifestFile = new File(runDir, "magma_geneset.manifest.json");
        writeGeneSetManifest(manifestFile, result.sets.size());
        result.manifestPath = manifestFile.getAbsolutePath();

        return result;
    }

    // ── Parsers (package-visible/static so tests can exercise them directly) ──────────

    /**
     * Parses MAGMA's .genes.out: a header line followed by whitespace-delimited rows
     * GENE CHR START STOP NSNPS NPARAM N ZSTAT P. Lines starting with '#' (MAGMA sometimes
     * writes a settings/comment preamble) are skipped, as are blank lines.
     */
    static List<GeneResult> parseGenesOut(File f) throws IOException {
        List<GeneResult> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            boolean sawHeader = false;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                if (!sawHeader) { sawHeader = true; continue; } // header row: GENE CHR START STOP ...
                String[] cols = line.split("\\s+");
                if (cols.length < 9) continue;
                GeneResult g = new GeneResult();
                g.gene = cols[0];
                g.chr = cols[1];
                g.start = Long.parseLong(cols[2]);
                g.stop = Long.parseLong(cols[3]);
                g.nsnps = Integer.parseInt(cols[4]);
                g.nparam = Integer.parseInt(cols[5]);
                g.n = Integer.parseInt(cols[6]);
                g.zstat = Double.parseDouble(cols[7]);
                g.p = Double.parseDouble(cols[8]);
                out.add(g);
            }
        }
        return out;
    }

    /**
     * Parses MAGMA's .gsa.out: a header line followed by whitespace-delimited rows
     * VARIABLE TYPE NGENES BETA BETA_STD SE P (competitive gene-set test).
     */
    static List<GeneSetResult> parseGsaOut(File f) throws IOException {
        List<GeneSetResult> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            boolean sawHeader = false;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                if (!sawHeader) { sawHeader = true; continue; } // header row: VARIABLE TYPE NGENES ...
                String[] cols = line.split("\\s+");
                if (cols.length < 7) continue;
                GeneSetResult r = new GeneSetResult();
                r.setId = cols[0];
                r.nGenes = Integer.parseInt(cols[2]);
                r.beta = Double.parseDouble(cols[3]);
                r.betaStd = Double.parseDouble(cols[4]);
                r.se = Double.parseDouble(cols[5]);
                r.p = Double.parseDouble(cols[6]);
                out.add(r);
            }
        }
        return out;
    }

    // ── Manifests ───────────────────────────────────────────────────────────────────────

    private static void writeGeneManifest(File f, int nGenes, int snpsWritten, int droppedNoRefId,
                                           int droppedStats, boolean perSnpN, int sampleN) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
            pw.print("{\"schema_version\":\"1.0\",\"method\":\"magma_gene_analysis\",\"method_version\":\"1.0\",");
            pw.printf("\"parameters\":{\"n_source\":\"%s\",\"sample_size\":%d},",
                perSnpN ? "per_snp" : "global", sampleN);
            pw.print("\"columns\":[");
            pw.print(String.join(",",
                colDef("gene", "string"), colDef("chr", "string"), colDef("start", "int"),
                colDef("stop", "int"), colDef("nsnps", "int"), colDef("nparam", "int"),
                colDef("n", "int"), colDef("zstat", "double"), colDef("p", "double")));
            pw.print("],");
            pw.printf("\"n_genes\":%d,\"snps_written\":%d,\"snps_dropped_no_ref_id\":%d,\"snps_dropped_stats\":%d,",
                nGenes, snpsWritten, droppedNoRefId, droppedStats);
            pw.printf("\"created_at\":%d}", System.currentTimeMillis());
        }
    }

    private static void writeGeneSetManifest(File f, int nSets) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
            pw.print("{\"schema_version\":\"1.0\",\"method\":\"magma_geneset_analysis\",\"method_version\":\"1.0\",");
            pw.print("\"parameters\":{},");
            pw.print("\"columns\":[");
            pw.print(String.join(",",
                colDef2("set_id", "string"), colDef2("n_genes", "int"), colDef2("beta", "double"),
                colDef2("beta_std", "double"), colDef2("se", "double"), colDef2("p", "double")));
            pw.print("],");
            pw.printf("\"n_sets\":%d,\"created_at\":%d}", nSets, System.currentTimeMillis());
        }
    }

    private static String colDef(String name, String type) {
        return String.format(
            "{\"name\":\"%s\",\"type\":\"%s\",\"scope\":\"per_gene\",\"method\":\"magma_gene_analysis\",\"method_version\":\"1.0\"}",
            name, type);
    }

    private static String colDef2(String name, String type) {
        return String.format(
            "{\"name\":\"%s\",\"type\":\"%s\",\"scope\":\"per_gene_set\",\"method\":\"magma_geneset_analysis\",\"method_version\":\"1.0\"}",
            name, type);
    }

    // ── Small helpers ───────────────────────────────────────────────────────────────────

    private static int findColumn(String header, String name) {
        if (header == null) return -1;
        String[] cols = header.split("\t", -1);
        for (int i = 0; i < cols.length; i++) {
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        }
        return -1;
    }

    private static int runProcess(List<String> cmd, File workDir, File logFile) throws IOException {
        System.out.printf("[MagmaAdapter] Running: %s%n", String.join(" ", cmd));
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd).directory(workDir).redirectErrorStream(true);
            Process proc = pb.start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                 PrintWriter log = new PrintWriter(new BufferedWriter(new FileWriter(logFile)))) {
                String line;
                while ((line = br.readLine()) != null) log.println(line);
            }
            return proc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("MAGMA process interrupted", e);
        }
    }

    private static String tailOf(File logFile, int maxLines) {
        if (!logFile.exists()) return null;
        try {
            List<String> lines = Files.readAllLines(logFile.toPath());
            int from = Math.max(0, lines.size() - maxLines);
            return String.join("\n", lines.subList(from, lines.size()));
        } catch (IOException e) {
            return null;
        }
    }
}
