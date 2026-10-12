# YouTube upload

**File:** `launch/video/lynxgwas-intro.mp4` (re-render: `python launch/video/make_video.py`)
**Thumbnail:** `launch/video/youtube-thumbnail.png` (1280×720)

## Title (≤ 100)
LYNXgwas: fine-mapping, heritability and causal-SNP ranking for your GWAS, on your own computer

## Description
One GWAS gives you hundreds of significant loci. Which variant is causal? Which gene does it act through? Does it replicate? Answering that usually takes a dozen tools, a dozen file formats and weeks of scripts, or uploading your data to a web server.

LYNXgwas does it in one app that runs on your computer. Point it at your summary statistics and it finds the loci and runs every analysis:
🔍 An interactive locus viewer: genome-wide skyline, regional plot, genes, regulatory tracks and an LD triangle on the same axis
🎯 Fine-mapping at every locus: SuSiE, Wakefield ABF and GCTA-COJO, with credible sets on the plot
🔗 Joint fine-mapping across datasets (mvSuSiE)
🏆 A causal-SNP ranking that combines all the evidence, with its reasons
🧬 SNP-heritability (built-in LDSC), local heritability (HESS), QC and a sample-size check
🧪 MAGMA at every locus and regulatory-element enrichment
🧫 Enhancer–gene links (scE2G) for 40 immune and pancreatic-islet cell types, with a per-cell-type enrichment test
🗺️ Many GWAS at once: the Locus Serpent Plot, the Gene Constellation and a search across every dataset
🤖 A built-in AI agent, with a local model by default
🔒 Private by design: nothing is uploaded

Robust and fast: a 7.66-million-SNP GWAS is read in 4 seconds and all 177 of its loci are processed in 37 minutes on one desktop; the built-in LDSC matches the reference ldsc to within a tenth of a standard error; wrong sample sizes are caught from the standard errors (9 of the 30 datasets).

It matches the published results: schizophrenia SNP-heritability 0.22 (published 0.24), type 2 diabetes 0.204 (published 0.196), and known causal variants recovered (PCSK9 R46L, SORT1 rs12740374, SH2B3 R262W). And it finds what only a multi-dataset tool can: the ranking's top SNP is the causal one at 60.8% of loci (50.3% for SuSiE alone), and joint fine-mapping narrows the credible set in 47% of cases.

Every example in the video comes from the article's 30 public GWAS across six diseases, which install into LYNXgwas in one click. No install? LYNXgwas Server runs the same tool in the browser, with public datasets and accounts for your own projects.

Free and open source (MIT).

⏱ Chapters
0:00 One GWAS, 7.6 million SNPs
0:13 The problem
0:22 Meet LYNXgwas
0:29 Three easy steps
0:47 The locus viewer
1:06 Fine-mapping: SuSiE, ABF, COJO, mvSuSiE
1:25 Causal-SNP ranking
1:45 Heritability and QC
1:59 Genes, enhancer–gene links and regulatory enrichment
2:23 Across datasets: Serpent Plot, Gene Constellation, search
2:48 The AI agent
2:55 Robust and fast
3:11 Matches the published results
3:26 What it found
3:42 The example data: 30 GWAS, 6 diseases
3:52 Use it online: LYNXgwas Server
4:07 Get it

🔗 Links
GitHub: https://github.com/AlsammanAlsamman/LYNXgwas
LYNXgwas Server: https://github.com/AlsammanAlsamman/lynxgwas-server
PyPI: https://pypi.org/project/lynxgwas/
Example data (Zenodo): <add the DOI once the record is published>

Install: pip install lynxgwas   →   lynxgwas

#GWAS #genetics #bioinformatics #finemapping #opensource

## Tags
LYNXgwas, GWAS, genome-wide association study, fine-mapping, SuSiE, mvSuSiE, MAGMA, LD score regression, heritability, causal variant, PCSK9, bioinformatics, genetics, statistical genetics, open source, python

## Settings
- Category: Science & Technology
- Audience: not made for kids
- Music: original, synthesized for this video (no copyright claims)

## After uploading: the README section
Paste this into `README.md` after the paragraph about LYNXgwas Server, with the video's id in place of
`YOUTUBE_VIDEO_ID` (the part after `watch?v=` in its link). The image is `docs/images/video-thumbnail.png`
(made by `make_thumbnail.py`, with a play button).

```markdown
## Watch the introduction

<p align="center">
  <a href="https://www.youtube.com/watch?v=YOUTUBE_VIDEO_ID">
    <img src="https://raw.githubusercontent.com/AlsammanAlsamman/LYNXgwas/main/docs/images/video-thumbnail.png"
         alt="LYNXgwas introduction video (4 minutes)" width="720">
  </a>
</p>

Four minutes on what LYNXgwas does, shown on the article's 30 GWAS datasets: the locus viewer, fine-mapping,
the causal-SNP ranking, heritability and QC, cross-dataset views, speed, agreement with published results,
and the online LYNXgwas Server.
```
