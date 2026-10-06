# YouTube upload

**File:** `launch/video/lynxgwas-intro.mp4` (re-render: `python launch/video/make_video.py`)
**Thumbnail:** `launch/video/youtube-thumbnail.png` (1280×720)

## Title (≤ 100)
LYNXgwas: from GWAS summary statistics to the causal variant, on your own computer (30 GWAS, 6 diseases)

## Description
One GWAS gives you hundreds of significant loci. Which variant is causal? Which gene does it act through? Does it replicate in other studies? Answering that usually takes a dozen tools, a dozen file formats and weeks of scripts, or uploading your data to a web server.

LYNXgwas does it all in one app that runs on your computer:
🎯 Fine-mapping at every locus: SuSiE, Wakefield ABF and GCTA-COJO
🔗 Joint fine-mapping across datasets (mvSuSiE)
🏆 A causal-SNP ranking that combines all the evidence
🧬 MAGMA per locus, LD-score and local (HESS) heritability
🗺️ Across datasets: the Locus Serpent Plot and the Gene Constellation
🔒 Private by design: nothing is uploaded

Shown here on the article's corpus of 30 public GWAS across six diseases (schizophrenia, type 2 diabetes, coronary artery disease, Alzheimer's disease, rheumatoid arthritis and inflammatory bowel disease; 2,387 loci). At PCSK9, every method agrees on rs11591147 (R46L) and the ranking scores it 0.993. Across the corpus the model's top SNP is the causal one at 60.8% of loci, against 50.3% for SuSiE alone.

Also online: LYNXgwas Server runs the same app in the browser for a lab or the public, with read-only public datasets and accounts for your own projects. The article's data can be added in one click from Zenodo, and projects move between computers as one .zip.

Free and open source (MIT).

⏱ Chapters
0:00 One GWAS, 7.6 million SNPs
0:15 The problem
0:26 Meet LYNXgwas
0:33 The test: 30 GWAS, 6 diseases
0:48 What it does
1:04 All 30 datasets
1:12 A locus close up: PCSK9
1:34 TCF7L2, IL23R, DRD2, SORT1
1:48 Locus Serpent Plot
2:04 Gene Constellation
2:13 Every analysis at a glance
2:25 The results
2:39 Online: LYNXgwas Server
2:54 Example data, Export and Import
3:07 Get it

🔗 Links
GitHub: https://github.com/AlsammanAlsamman/LYNXgwas
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
         alt="LYNXgwas introduction video (3 minutes)" width="720">
  </a>
</p>

Three minutes on the article's 30 GWAS datasets: fine-mapping, joint fine-mapping, the causal-SNP ranking,
the Locus Serpent Plot and Gene Constellation, the online LYNXgwas Server, and the one-click example data.
```
