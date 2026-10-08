# Enhancer–gene links (scE2G)

LYNXgwas ships with the enhancer–gene links of **scE2G** for 40 cell types: 9 from peripheral blood (PBMC),
22 from bone marrow and 9 from pancreatic islets. In every locus they show which promoters and enhancers lie
near the lead SNP, which genes those elements regulate, and in which cell types; the project summary page
tests whether the lead SNPs sit in these elements more often than random SNPs, per cell type.

## Source and citation

scE2G is an enhancer–gene prediction model trained on CRISPR perturbation data and applied to single-cell
multiome (RNA + ATAC) data. Please cite it when you use these links:

> Sheth MU, Qiu W-L, et al. Mapping enhancer-gene regulatory interactions from single-cell data.
> bioRxiv 2024.11.23.624931 (Engreitz lab). Predictions: https://data.broadinstitute.org/scE2Gpredictions/

LYNXgwas uses the public `scE2G_multiome` predictions of the `PBMC_9_Cell_Types`, `BMMC_22_Cell_Types` and
`Pancreatic_Islets` folders (the "merged" files are not used).

## What is bundled, and how it was reduced

Each public file lists every candidate element–gene pair genome-wide (about 10 million pairs, 27 columns,
~700 MB per cell type, ~29 GB for the 40). LYNXgwas keeps, per cell type:

- only the pairs scE2G itself calls: **E2G.Score.qnorm ≥ 0.177**, the scE2G_multiome threshold
  (about 60,000–75,000 pairs per cell type);
- the element (coordinates and class: promoter, genic or intergenic enhancer), the target gene and its TSS,
  whether the link is the gene's own promoter, and the scores E2G (normalized and raw), ABC, Kendall and
  ARC-E2G, rounded to 4 significant digits;
- both genome builds: the original GRCh38, and GRCh37 lifted with the UCSC `hg38ToHg19` chain (pairs whose
  element or TSS does not lift cleanly, about 0.3%, are left out of both). GRCh38 is stored as offsets from
  GRCh37 (`d38`, `dtss38`, `dlen38`: start38 = start + d38, end38 = start38 + (end − start) + dlen38,
  tss38 = tss + d38 + dtss38), which is exact and about a quarter smaller than absolute positions.

The model's input features (ATAC, RNA, enhancer counts) are dropped. The result is about 1.7 MB per cell type,
69 MB in all (2.73 million links), in `regulatory_links/sce2g/` with a `manifest.json` (labels, tissue,
lineage, counts, checksums and the source URL of each file). `scripts/regulatory/build_sce2g_pack.py` rebuilds
it from the public files (one file at a time, deleting each after compaction), then
`scripts/regulatory/delta38.py regulatory_links/sce2g` converts the GRCh38 columns to offsets;
`scripts/regulatory/compact_sce2g.py` compacts one file.

## In the locus viewer

A track under the histone-mark tracks shows:

- the elements of the selected cell types, merged across cell types into regions; colour = class
  (promoter, genic enhancer, intergenic enhancer), bar height = number of cell types calling the element;
- the element nearest the lead SNP, outlined, with its distance;
- arcs from that element to the TSSs of its target genes, width = E2G score;
- credible-set SNPs (latest SuSiE run) that sit inside an element.

Under the plot, a panel gives the cell-type selection (immune cells, PBMC, bone marrow, pancreatic islets or
all), the nearest element and its target genes with E2G, ABC and Kendall scores, the genes linked to
credible-set SNPs ranked by Σ PIP × E2G, and the nearest element in every cell type. The default selection
follows the project's disease: pancreatic islets for type 2 diabetes, immune cells for immune diseases, all
cell types otherwise.

## On the summary page

- **Where the lead SNPs sit:** inside an element, within 1 kb, within 10 kb or farther; the class of the
  nearest element; the cumulative distance against matched random SNPs.
- **Enrichment:** per cell type and for any selected cell type, the number of lead SNPs inside / within 1 kb
  / within 10 kb against 1,000 sets of random SNPs matched on the signed distance to the nearest gene TSS
  (each lead SNP's offset placed at a random TSS of the project's gene annotation). One-sided empirical
  p = (1 + sets ≥ observed) / (1 + sets); q by Benjamini–Hochberg across cell types. Heatmap, sortable table
  and CSV download.
- **Credible sets:** loci with credible-set SNPs inside elements and the genes they point to.
- **Every lead SNP:** distance, class, cell types with the SNP inside, linked gene.

MHC loci are excluded. Results are cached in the project folder per cell-type selection.

## Endpoints

`GET /api/sce2g/info`, `GET /api/project/{id}/locus/{n}/sce2g?cts=<preset or ids>`,
`GET /api/project/{id}/sce2g-stats?cts=<preset or ids>[&refresh=1]`.
