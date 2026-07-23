package export;

/**
 * Minimal copy of the project's column mapping + raw GWAS file path, passed in
 * from the root-package Config (export/ can't import the default-package Config
 * class directly) so the exporter can re-read the raw GWAS file for columns
 * that aren't part of the standardized pipeline (e.g. METAL's Direction/HetISq/
 * HetChiSq/HetDf/HetPVal, MR-MEGA's per-cohort betas/heterogeneity stats).
 */
public class GwasSchema {
    public String gwasFile;
    public String colChr, colPos, colEa, colNea, colPvalue;
    public String colRsid, colVarid, colBeta, colOr, colSe, colN, colMaf, colInfo;
}
