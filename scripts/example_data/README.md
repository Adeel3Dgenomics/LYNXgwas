# LYNXgwas article data: 30 GWAS datasets, 6 diseases (results only)

This is the processed corpus behind the LYNXgwas article: 30 public GWAS summary-statistics datasets, five
per disease for schizophrenia, type 2 diabetes, coronary artery disease, Alzheimer's disease, rheumatoid
arthritis and inflammatory bowel disease, all on GRCh37, with 2,387 genome-wide-significant loci.

## Using it

The easiest way is from LYNXgwas itself: on the home page, click **Add example data** (or **Download** on
the banner shown when the example data is not installed yet). LYNXgwas downloads this file, checks it,
merges the 30 projects and the two cross-dataset runs into your workspace, and fetches the GENCODE gene
annotation the projects use. Your own projects are never touched.

To install a copy you downloaded yourself, use **Add example data → Install from a file**.

## What is included

For each of the 30 datasets, everything LYNXgwas shows:

- the loci, the genome-wide and per-locus viewer data, including the LD shown in the viewer
- fine-mapping (SuSiE, Wakefield ABF, GCTA-COJO), the MAGMA gene test, local heritability (HESS),
  LD-score heritability, regulatory-element enrichment, the sample-size check and the causal-SNP ranking,
  as the latest run of each tool at each locus, with its provenance
- the two cross-dataset runs over all 30 datasets (the Locus Serpent Plot and the Gene Constellation), with
  their joint fine-mapping (mvSuSiE) runs
- the Roadmap Epigenomics histone peaks used for the regulatory tracks

## What is not included, and what that means

- **The GWAS summary statistics.** They belong to the consortia that published them (see the table below)
  and are not redistributed here. Every result can be browsed without them. Re-running an analysis or
  defining new loci needs the original file: download it from its source and point the project at it.
- **The LD matrices, reference-panel subsets and other intermediate files** that LYNXgwas rebuilds itself
  when an analysis is run again. They need the 1000 Genomes reference panel, which is not included either.

The projects are marked *example data (results only)* in LYNXgwas, and the actions that need the original
summary statistics are switched off for them.

## Datasets

| Disease | Dataset | LYNXgwas project | Loci |
|---|---|---|---|
| Schizophrenia | Trubetskoy et al. 2022 (PGC3) | scz-demo | 177 |
| | Pardiñas et al. 2018 (CLOZUK) | scz-src-clozuk2018 | 140 |
| | PGC 2011 (PGC1) | scz-src-pgc1-2011 | 8 |
| | Ripke et al. 2013 (Sweden) | scz-src-sweden2013 | 10 |
| | Ripke et al. 2014 (PGC2) | scz-src-pgc2-2014 | 94 |
| Type 2 diabetes | Xue et al. 2018 | t2d-demo | 135 |
| | DIAMANTE 2022 (Mahajan et al.) | t2d-src-diamante2022 | 212 |
| | FinnGen R12 | t2d-src-finngen-r12 | 322 |
| | Mahajan et al. 2018 | t2d-src-mahajan2018 | 58 |
| | Scott et al. 2017 (DIAGRAM) | t2d-src-scott2017 | 41 |
| Coronary artery disease | van der Harst & Verweij 2018 | cad-demo | 75 |
| | Aragam et al. 2022 | cad-src-aragam2022 | 209 |
| | Hartiala et al. 2021 | cad-src-hartiala2021 | 81 |
| | Nam et al. 2026 | cad-src-nam2026 | 72 |
| | Nikpay et al. 2015 | cad-src-nikpay2015 | 42 |
| Alzheimer's disease | Kunkle et al. 2019 | alz-demo | 15 |
| | Bellenguez et al. 2022 | alz-src-bellenguez2022 | 57 |
| | EADB 2026 | alz-src-eadb2026 | 65 |
| | FinnGen R12 | alz-src-finngen-r12 | 19 |
| | Nicolas et al. 2025 | alz-src-nicolas2025 | 44 |
| Rheumatoid arthritis | Ishigaki et al. 2022 | ra-demo | 65 |
| | Eyre et al. 2012 | ra-src-eyre2012 | 13 |
| | FinnGen R12 | ra-src-finngen-r12 | 29 |
| | Glanville et al. 2021 | ra-src-glanville2021 | 10 |
| | Stahl et al. 2010 | ra-src-stahl2010 | 5 |
| Inflammatory bowel disease | de Lange et al. 2017 (Crohn's disease) | ibd-demo | 91 |
| | de Lange et al. 2017 (ulcerative colitis) | ibd-src-delange2017uc | 64 |
| | FinnGen R12 | ibd-src-finngen-r12 | 89 |
| | Liu et al. 2015 | ibd-src-liu2015 | 138 |
| | Zorina-Lichtenwalter et al. 2023 | ibd-src-zorina2023 | 7 |

If you use these results, please cite the LYNXgwas article and the original publication of each dataset
you use.

## Files

- `projects/<project>/` — one LYNXgwas project per dataset
- `output/multi_locus/<run>/` and `output/mvsusie/<run>/` — the cross-dataset runs and their joint
  fine-mapping
- `regulatory_data/` — Roadmap Epigenomics histone peaks (narrowPeak)
- `pack.json` — the pack's contents, and the files whose paths LYNXgwas fills in when it installs them

## Licence

The processed results are released under CC BY 4.0 (see `LICENSE.txt`). The underlying GWAS summary
statistics, the Roadmap Epigenomics data and GENCODE remain under their own terms.
