import java.io.*;
import java.util.*;
import loci.LociIdentifier;

/**
 * Locus Matrix phase 3: for each selected dataset, streams its RAW GWAS file
 * once (not the merged file) to compute per-locus running aggregates: the
 * best (lowest) p SNP with its beta/OR/se, and counts of SNPs with
 * p&lt;5e-5 and p&lt;5e-8. Deliberately does NOT reuse GwasParser.parse(),
 * whose downsampling (keep p&lt;1e-4, random-sample the rest) can drop the
 * true best SNP for a dataset with no strong signal in a large/dense locus.
 * This sweep only tracks scalar aggregates, so nothing needs downsampling.
 */
public class MultiLocusScanner {

    private static class LocusAgg {
        String bestSnpId = "";
        String bestChr   = "";
        long   bestPos   = 0;
        double bestP     = Double.NaN;
        String ea        = "";
        String nea       = "";
        double beta      = Double.NaN;
        double or        = Double.NaN;
        double se        = Double.NaN;
        double maf       = Double.NaN;
        double n         = Double.NaN;
        double info      = Double.NaN;
        int    nP5e5     = 0;
        int    nP5e8     = 0;
        Map<String, String> extra = Collections.emptyMap();
    }

    /** Result of scanning one dataset's file: per-locus aggregates plus the list of that
     *  dataset's own header columns not already mapped to a structured field (e.g. METAL's
     *  Direction/HetISq/HetChiSq/HetDf/HetPVal, or MR-MEGA's per-cohort betas/het stats). */
    private static class ScanOneResult {
        Map<Integer, LocusAgg> agg;
        List<String> extraColumns;
    }

    public static MultiLocusResult scan(List<Config> configs, List<String> datasetIds,
                                         List<LociIdentifier.IdentifiedLocus> loci,
                                         GffParser gff, MultiLocusProgress progress) throws IOException {
        progress.phase = "scan";
        progress.datasetTotal = configs.size();

        List<LociIdentifier.IdentifiedLocus> sorted = new ArrayList<>(loci);
        sorted.sort(Comparator.comparingInt((LociIdentifier.IdentifiedLocus l) -> Locus.chrToInt(l.chr))
                              .thenComparingLong(l -> l.start));

        Map<String, Map<Integer, LocusAgg>> perDataset = new LinkedHashMap<>();
        Map<String, List<String>> extraColumnsByDataset = new LinkedHashMap<>();
        for (int di = 0; di < configs.size(); di++) {
            Config cfg = configs.get(di);
            progress.datasetIndex = di;
            progress.currentDataset = datasetIds.get(di);
            ScanOneResult sr = scanOne(cfg, sorted);
            perDataset.put(datasetIds.get(di), sr.agg);
            extraColumnsByDataset.put(datasetIds.get(di), sr.extraColumns);
        }
        progress.datasetIndex = configs.size();

        MultiLocusResult result = new MultiLocusResult();
        for (int di = 0; di < configs.size(); di++) {
            Config cfg = configs.get(di);
            String effectType = !cfg.colOr.isEmpty() ? "OR" : "Beta";
            MultiLocusResult.DatasetInfo di_ = new MultiLocusResult.DatasetInfo(
                datasetIds.get(di), new File(cfg.outputDir).getName(), effectType);
            di_.extraColumns = extraColumnsByDataset.get(datasetIds.get(di));
            result.datasets.add(di_);
        }

        for (int li = 0; li < sorted.size(); li++) {
            LociIdentifier.IdentifiedLocus loc = sorted.get(li);
            MultiLocusResult.LocusRow row = new MultiLocusResult.LocusRow();
            row.index  = li;
            row.chr    = loc.chr;
            row.start  = loc.start;
            row.end    = loc.end;
            row.sizeBp = loc.end - loc.start + 1;

            String bestDs = null;
            double bestP = Double.POSITIVE_INFINITY;
            String bestSnp = "";
            long bestPos = 0;
            String bestChr = loc.chr;

            for (String dsId : datasetIds) {
                LocusAgg a = perDataset.get(dsId).get(li);
                MultiLocusResult.DatasetLocusStat stat = new MultiLocusResult.DatasetLocusStat();
                if (a != null && !Double.isNaN(a.bestP)) {
                    stat.bestSnpId = a.bestSnpId;
                    stat.bestChr   = a.bestChr;
                    stat.bestPos   = a.bestPos;
                    stat.bestP     = a.bestP;
                    stat.ea        = a.ea;
                    stat.nea       = a.nea;
                    stat.beta      = a.beta;
                    stat.or        = a.or;
                    stat.se        = a.se;
                    stat.maf       = a.maf;
                    stat.n         = a.n;
                    stat.info      = a.info;
                    stat.nP5e5     = a.nP5e5;
                    stat.nP5e8     = a.nP5e8;
                    stat.extra     = a.extra;
                    if (a.bestP < bestP) {
                        bestP = a.bestP;
                        bestDs = dsId;
                        bestSnp = a.bestSnpId;
                        bestPos = a.bestPos;
                        bestChr = a.bestChr;
                    }
                }
                row.cells.put(dsId, stat);
            }

            row.overallBestDatasetId = bestDs != null ? bestDs : "";
            row.overallBestSnpId     = bestSnp;
            row.overallBestP         = bestDs != null ? bestP : Double.NaN;

            if (bestDs != null && gff != null) {
                List<Gene> genes = new ArrayList<>(gff.overlapping(loc.chr, loc.start, loc.end));
                List<String> names = gff.nearestGeneNames(bestPos, genes); // sorts `genes` by distance in place
                if (!names.isEmpty()) {
                    row.nearestGene = String.join(", ", names);
                    row.geneDist = genes.get(0).distanceTo(bestPos);
                }
            }

            result.loci.add(row);
            progress.lociFound = result.loci.size();
        }

        return result;
    }

    private static ScanOneResult scanOne(Config cfg,
            List<LociIdentifier.IdentifiedLocus> sorted) throws IOException {
        Map<Integer, LocusAgg> agg = new HashMap<>();
        ScanOneResult result = new ScanOneResult();
        result.agg = agg;
        result.extraColumns = new ArrayList<>();
        List<Integer> activeIdx = new ArrayList<>();
        int nextLocusIdx = 0;
        String currentChr = null;

        try (BufferedReader br = new BufferedReader(new FileReader(cfg.gwasFile), 1024 * 1024)) {
            String header = br.readLine();
            if (header == null) return result;
            String[] cols = header.trim().split("\t");
            int iChr   = GwasParser.colIdx(cols, cfg.colChr);
            int iPos   = GwasParser.colIdx(cols, cfg.colPos);
            int iP     = GwasParser.colIdx(cols, cfg.colPvalue);
            int iRsid  = cfg.colRsid.isEmpty()  ? -1 : GwasParser.colIdx(cols, cfg.colRsid);
            int iVarid = cfg.colVarid.isEmpty() ? -1 : GwasParser.colIdx(cols, cfg.colVarid);
            int iEa    = cfg.colEa.isEmpty()    ? -1 : GwasParser.colIdx(cols, cfg.colEa);
            int iNea   = cfg.colNea.isEmpty()   ? -1 : GwasParser.colIdx(cols, cfg.colNea);
            int iBeta  = cfg.colBeta.isEmpty()  ? -1 : GwasParser.colIdx(cols, cfg.colBeta);
            int iOr    = cfg.colOr.isEmpty()    ? -1 : GwasParser.colIdx(cols, cfg.colOr);
            int iSe    = cfg.colSe.isEmpty()    ? -1 : GwasParser.colIdx(cols, cfg.colSe);
            int iMaf   = cfg.colMaf.isEmpty()   ? -1 : GwasParser.colIdx(cols, cfg.colMaf);
            int iN     = cfg.colN.isEmpty()     ? -1 : GwasParser.colIdx(cols, cfg.colN);
            int iInfo  = cfg.colInfo.isEmpty()  ? -1 : GwasParser.colIdx(cols, cfg.colInfo);
            if (iChr < 0 || iPos < 0 || iP < 0)
                throw new IOException("Missing required columns in " + cfg.gwasFile);

            // Every header column not already mapped to a structured field becomes a
            // generic "extra" column (e.g. METAL's Direction/HetISq/HetChiSq/HetDf/HetPVal,
            // MR-MEGA's per-cohort betas/heterogeneity stats, etc.), captured for the best SNP.
            Set<Integer> mappedIdx = new HashSet<>(Arrays.asList(
                iChr, iPos, iP, iRsid, iVarid, iEa, iNea, iBeta, iOr, iSe, iMaf, iN, iInfo));
            List<Integer> extraIdx = new ArrayList<>();
            for (int i = 0; i < cols.length; i++) {
                if (!mappedIdx.contains(i)) {
                    extraIdx.add(i);
                    result.extraColumns.add(cols[i].trim());
                }
            }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = GwasParser.splitTab(line);
                if (f.length <= Math.max(iChr, Math.max(iPos, iP))) continue;

                String chr = f[iChr].trim();
                long pos; double p;
                try {
                    pos = Long.parseLong(f[iPos].trim());
                    p   = Double.parseDouble(f[iP].trim());
                } catch (NumberFormatException e) { continue; }

                int chrI = Locus.chrToInt(chr);
                if (!chr.equals(currentChr)) {
                    activeIdx.clear();
                    currentChr = chr;
                    while (nextLocusIdx < sorted.size()
                           && Locus.chrToInt(sorted.get(nextLocusIdx).chr) < chrI) {
                        nextLocusIdx++;
                    }
                }

                while (nextLocusIdx < sorted.size()) {
                    LociIdentifier.IdentifiedLocus cand = sorted.get(nextLocusIdx);
                    if (Locus.chrToInt(cand.chr) != chrI) break;
                    if (cand.start > pos) break;
                    activeIdx.add(nextLocusIdx);
                    nextLocusIdx++;
                }
                activeIdx.removeIf(idx -> sorted.get(idx).end < pos);
                if (activeIdx.isEmpty()) continue;

                String rsid  = (iRsid  >= 0 && iRsid  < f.length) ? f[iRsid].trim()  : "";
                String varid = (iVarid >= 0 && iVarid < f.length) ? f[iVarid].trim() : chr + ":" + pos;
                String id    = (!rsid.isEmpty() && !rsid.equals(".")) ? rsid : varid;
                String ea    = (iEa   >= 0 && iEa   < f.length) ? f[iEa].trim()  : "";
                String nea   = (iNea  >= 0 && iNea  < f.length) ? f[iNea].trim() : "";
                double beta  = (iBeta >= 0 && iBeta < f.length) ? parseD(f[iBeta]) : Double.NaN;
                double or    = (iOr   >= 0 && iOr   < f.length) ? parseD(f[iOr])   : Double.NaN;
                double se    = (iSe   >= 0 && iSe   < f.length) ? parseD(f[iSe])   : Double.NaN;
                double maf   = (iMaf  >= 0 && iMaf  < f.length) ? parseD(f[iMaf])  : Double.NaN;
                double n     = (iN    >= 0 && iN    < f.length) ? parseD(f[iN])    : Double.NaN;
                double info  = (iInfo >= 0 && iInfo < f.length) ? parseD(f[iInfo]) : Double.NaN;

                for (int idx : activeIdx) {
                    LocusAgg a = agg.computeIfAbsent(idx, k -> new LocusAgg());
                    if (p < 5e-5) a.nP5e5++;
                    if (p < 5e-8) a.nP5e8++;
                    if (Double.isNaN(a.bestP) || p < a.bestP) {
                        a.bestP     = p;
                        a.bestSnpId = id;
                        a.bestChr   = chr;
                        a.bestPos   = pos;
                        a.ea        = ea;
                        a.nea       = nea;
                        a.beta      = beta;
                        a.or        = or;
                        a.se        = se;
                        a.maf       = maf;
                        a.n         = n;
                        a.info      = info;
                        Map<String, String> extra = new LinkedHashMap<>();
                        for (int ci = 0; ci < extraIdx.size(); ci++) {
                            int idx2 = extraIdx.get(ci);
                            extra.put(result.extraColumns.get(ci), idx2 < f.length ? f[idx2].trim() : "");
                        }
                        a.extra = extra;
                    }
                }
            }
        }
        return result;
    }

    private static double parseD(String v) {
        v = v.trim();
        if (v.isEmpty() || v.equals(".") || v.equalsIgnoreCase("NA") || v.equalsIgnoreCase("nan"))
            return Double.NaN;
        try { return Double.parseDouble(v); }
        catch (NumberFormatException e) { return Double.NaN; }
    }
}
