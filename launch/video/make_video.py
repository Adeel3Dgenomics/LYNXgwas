"""Render the LYNXgwas intro video (1080p, about 4 minutes) with an original synthesized soundtrack.

    pip install pillow numpy imageio-ffmpeg
    python launch/video/capture.py && python launch/video/capture_extra.py && python launch/video/capture_features.py
    python launch/video/make_video.py            ->  launch/video/lynxgwas-intro.mp4
    python launch/video/make_video.py --preview  ->  a few still frames only (fast)

Every plot and number comes from the article's 30-dataset corpus: the Manhattan plots are drawn from the
projects' genome skylines (launch/video/data/skylines.json) and the screens are real LYNXgwas views.
"""

from __future__ import annotations

import json
import math
import random
import subprocess
import sys
import wave
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
SHOTS = HERE / "shots"
ICON = ROOT / "docs" / "images" / "icon.png"
SKY = json.loads((HERE / "data" / "skylines.json").read_text())
W, H, FPS = 1920, 1080, 30

BG = (7, 15, 19)
TEXT = (238, 243, 242)
MUTED = (160, 178, 178)
TEAL = (15, 110, 86)
MINT = (52, 211, 170)
AMBER = (250, 199, 117)
RED = (239, 68, 68)
ORANGE = (251, 146, 60)
SKYBLUE = (125, 211, 252)
VIOLET = (167, 139, 250)
CHR_A, CHR_B = (40, 140, 120), (70, 95, 160)


# --------------------------------------------------------------------------- helpers
def font(size: int, bold: bool = False, mono: bool = False) -> ImageFont.FreeTypeFont:
    name = "consola.ttf" if mono else ("segoeuib.ttf" if bold else "segoeui.ttf")
    return ImageFont.truetype(f"C:/Windows/Fonts/{name}", size)


EMOJI = ImageFont.truetype("C:/Windows/Fonts/seguiemj.ttf", 109)


def ease(x: float) -> float:
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def fade(t: float, start: float, end: float, f: float = 0.6) -> float:
    return ease((t - start) / f) * ease((end - t) / f)


_img: dict = {}


def img(name: str) -> Image.Image:
    if name not in _img:
        p = SHOTS / name if not Path(name).is_absolute() else Path(name)
        _img[name] = Image.open(p).convert("RGBA")
    return _img[name]


def background(t: float) -> Image.Image:
    """Near-black with two slowly drifting glows (teal, amber) and a faint grid of 'genome' dots."""
    base = Image.new("RGB", (W, H), BG)
    glow = Image.new("RGB", (W // 8, H // 8), (0, 0, 0))
    g = ImageDraw.Draw(glow)
    for cx, cy, r, col in ((0.8 + 0.05 * math.sin(t / 7), 0.18, 75, TEAL),
                           (0.12, 0.72 + 0.05 * math.cos(t / 9), 55, (120, 80, 20))):
        g.ellipse((cx * W / 8 - r, cy * H / 8 - r, cx * W / 8 + r, cy * H / 8 + r), fill=col)
    glow = glow.filter(ImageFilter.GaussianBlur(30)).resize((W, H), Image.BILINEAR)
    return Image.blend(base, glow, 0.30)


def text_center(im, txt, y, fnt, color, alpha=1.0, x=None):
    if alpha <= 0:
        return
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    w = d.textlength(txt, font=fnt)
    d.text((((W if x is None else 2 * x) - w) / 2, y), txt, font=fnt, fill=color + (int(255 * alpha),))
    im.alpha_composite(layer)


def text_at(im, txt, xy, fnt, color, alpha=1.0):
    if alpha <= 0:
        return
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(layer).text(xy, txt, font=fnt, fill=color + (int(255 * alpha),))
    im.alpha_composite(layer)


def wrap(txt, fnt, width):
    d = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    lines, cur = [], ""
    for wd in txt.split():
        test = (cur + " " + wd).strip()
        if d.textlength(test, font=fnt) > width and cur:
            lines.append(cur)
            cur = wd
        else:
            cur = test
    lines.append(cur)
    return lines


def gradient_text(im, txt, y, fnt, alpha=1.0, c1=MINT, c2=AMBER, x=None):
    if alpha <= 0:
        return
    d = ImageDraw.Draw(im)
    w = int(d.textlength(txt, font=fnt))
    h = fnt.size + 24
    mask = Image.new("L", (w, h), 0)
    ImageDraw.Draw(mask).text((0, 0), txt, font=fnt, fill=int(255 * alpha))
    grad = Image.new("RGBA", (w, h))
    gd = ImageDraw.Draw(grad)
    for xx in range(w):
        k = xx / max(w - 1, 1)
        gd.line([(xx, 0), (xx, h)], fill=tuple(int(c1[i] + (c2[i] - c1[i]) * k) for i in range(3)) + (255,))
    grad.putalpha(mask)
    im.alpha_composite(grad, (int(((W if x is None else 2 * x) - w) / 2), y))


def emoji(im, chars, xy, size, alpha=1.0):
    if alpha <= 0:
        return
    chars = chars.replace(chr(0xFE0F), "")
    tmp = Image.new("RGBA", (int(EMOJI.size * 1.3) * len(chars), int(EMOJI.size * 1.3)), (0, 0, 0, 0))
    ImageDraw.Draw(tmp).text((0, 0), chars, font=EMOJI, embedded_color=True)
    tmp = tmp.crop(tmp.getbbox())
    k = size / EMOJI.size
    tmp = tmp.resize((max(1, int(tmp.width * k)), max(1, int(tmp.height * k))), Image.LANCZOS)
    if alpha < 1:
        tmp.putalpha(tmp.getchannel("A").point(lambda a: int(a * alpha)))
    im.alpha_composite(tmp, (int(xy[0]), int(xy[1])))


def frame(im, src, crop, box, alpha=1.0, radius=18, shadow=True, glow=None):
    """Draws the `crop` (x0, y0, x1, y1) of screenshot `src` into `box` (x, y, w, h), keeping its aspect."""
    if alpha <= 0:
        return
    s = img(src)
    region = s.crop(tuple(int(v) for v in crop))
    bx, by, bw, bh = box
    k = min(bw / region.width, bh / region.height)
    pic = region.resize((max(1, int(region.width * k)), max(1, int(region.height * k))), Image.BILINEAR)
    x, y = int(bx + (bw - pic.width) / 2), int(by + (bh - pic.height) / 2)
    m = Image.new("L", pic.size, 0)
    ImageDraw.Draw(m).rounded_rectangle((0, 0, pic.width - 1, pic.height - 1), radius, fill=int(255 * alpha))
    pic.putalpha(m)
    if shadow:
        sh = Image.new("RGBA", (pic.width + 140, pic.height + 140), (0, 0, 0, 0))
        ImageDraw.Draw(sh).rounded_rectangle((70, 90, pic.width + 70, pic.height + 90), radius + 6,
                                             fill=(glow or (0, 0, 0)) + (int(170 * alpha),))
        im.alpha_composite(sh.filter(ImageFilter.GaussianBlur(32)), (x - 70, y - 70))
    im.alpha_composite(pic, (x, y))
    return x, y, pic.width, pic.height


def lerp_box(a, b, k):
    k = ease(k)
    return tuple(a[i] + (b[i] - a[i]) * k for i in range(4))


def caption(im, title, sub, alpha, y=900):
    text_center(im, title, y, font(54, bold=True), TEXT, alpha)
    if sub:
        text_center(im, sub, y + 76, font(32), MUTED, alpha)


def section(im, t, kicker, title, sub, dur=3.4):
    a = ease(t / 0.6) * ease((dur - t) / 0.5)
    text_center(im, kicker, 330, font(38, bold=True), MINT, a)
    gradient_text(im, title, 385, font(112, bold=True), a)
    text_center(im, sub, 560, font(38), MUTED, a)


def pill(im, txt, xy, fnt, fg, bg, alpha=1.0, pad=(22, 10)):
    if alpha <= 0:
        return 0
    d = ImageDraw.Draw(im)
    w = d.textlength(txt, font=fnt)
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    x, y = xy
    ld.rounded_rectangle((x, y, x + w + 2 * pad[0], y + fnt.size + 2 * pad[1] + 6), 30, fill=bg + (int(235 * alpha),))
    ld.text((x + pad[0], y + pad[1] - 2), txt, font=fnt, fill=fg + (int(255 * alpha),))
    im.alpha_composite(layer)
    return w + 2 * pad[0]


# --------------------------------------------------------------------------- Manhattan plots
def chrom_index(sky):
    n, L = len(sky["bins"]), sky["genome_length"]
    offs = [c["offset"] for c in sky["chromosomes"]]
    idx = []
    for i in range(n):
        pos = (i + 0.5) / n * L
        j = 0
        while j + 1 < len(offs) and offs[j + 1] <= pos:
            j += 1
        idx.append(j)
    return idx


_chr_cache: dict = {}


def manhattan(im, pid, box, progress=1.0, alpha=1.0, cap=40.0, glow=True, label=None, small=False):
    """A Manhattan plot from a project's genome skyline: 3,000 bins of the strongest -log10(p) in each."""
    if alpha <= 0:
        return
    sky = SKY[pid]
    if pid not in _chr_cache:
        _chr_cache[pid] = chrom_index(sky)
    ci = _chr_cache[pid]
    bins = sky["bins"]
    n = len(bins)
    bx, by, bw, bh = box
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    thr = by + bh - 7.3 / cap * bh
    shown = int(n * min(max(progress, 0), 1))
    hot = []
    for i in range(shown):
        v = min(bins[i], cap)
        x = bx + i / n * bw
        top = by + bh - v / cap * bh
        if bins[i] >= 7.3:
            col = RED if bins[i] >= 20 else ORANGE
            hot.append((x, top, col))
        else:
            col = CHR_A if ci[i] % 2 == 0 else CHR_B
        d.line([(x, by + bh), (x, top)], fill=col + (int(235 * alpha),), width=1 if small else 2)
    # genome-wide significance line
    for xx in range(int(bx), int(bx + bw), 18):
        d.line([(xx, thr), (xx + 9, thr)], fill=RED + (int(170 * alpha),), width=2)
    if glow and hot:
        g = Image.new("RGBA", im.size, (0, 0, 0, 0))
        gd = ImageDraw.Draw(g)
        for x, top, col in hot:
            gd.ellipse((x - 9, top - 9, x + 9, top + 9), fill=col + (int(120 * alpha),))
        im.alpha_composite(g.filter(ImageFilter.GaussianBlur(9)))
    im.alpha_composite(layer)
    if label:
        text_at(im, label, (bx, by - (40 if small else 56)), font(30 if small else 40, bold=True), TEXT, alpha)


# --------------------------------------------------------------------------- scenes
def s_hook(t):
    """One GWAS draws itself: 7.66 million SNPs, and then the questions."""
    im = background(t).convert("RGBA")
    box = (110, 260, W - 220, 520)
    manhattan(im, "scz-demo", box, progress=ease(t / 6.0), alpha=ease(t / 0.8) * (1 - 0.55 * ease((t - 8.5) / 1.2)))
    n = int(SKY["scz-demo"]["snps"] * ease(t / 6.0))
    text_center(im, f"{n:,} SNPs", 120, font(80, bold=True, mono=True), TEXT, fade(t, 0.3, 7.0))
    text_center(im, "one schizophrenia GWAS", 225, font(34), MUTED, fade(t, 0.8, 7.0))
    text_center(im, "177 genome-wide significant loci.", 120, font(76, bold=True), AMBER, fade(t, 7.2, 10.2))
    q = [("Which variant is causal?", 8.4), ("Which gene does it act through?", 9.4), ("Does it replicate in other studies?", 10.4)]
    for i, (line, at) in enumerate(q):
        text_center(im, line, 760 + i * 64, font(50, bold=True), TEXT if i < 2 else MINT, fade(t, at, 13.0))
    return im


TOOLS = ["PLINK", "SuSiE", "FINEMAP / ABF", "GCTA-COJO", "MAGMA", "LDSC", "HESS", "mvSuSiE", "coloc", "SuSiEx",
         "GWAMA", "liftOver", "LocusZoom", "R scripts", "bash", "Excel", "reference panels", "GENCODE",
         "Roadmap peaks", "HapMap3", "CADD", "1000 Genomes", "python", "awk"]


def s_problem(t):
    im = background(t + 15).convert("RGBA")
    rnd = random.Random(3)
    collapse = ease((t - 9.0) / 1.4)
    for i, name in enumerate(TOOLS):
        x0 = rnd.uniform(140, W - 360)
        y0 = rnd.choice((rnd.uniform(90, 330), rnd.uniform(700, H - 120)))
        vx, vy = rnd.uniform(-22, 22), rnd.uniform(-6, 6)
        x = x0 + vx * t + math.sin(t * 1.3 + i) * 18
        y = y0 + vy * t + math.cos(t * 1.1 + i) * 14
        x = x + (W / 2 - 80 - x) * collapse
        y = y + (H / 2 - y) * collapse
        a = ease((t - i * 0.1) / 0.5) * (1 - collapse) * 0.9
        size = rnd.choice((34, 40, 46, 52))
        col = rnd.choice((TEXT, MUTED, SKYBLUE, AMBER, MINT))
        text_at(im, name, (x, y), font(size, bold=size > 40, mono=rnd.random() < 0.3), col, a)
    band = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(band).rectangle((0, 400, W, 690), fill=BG + (int(240 * fade(t, 0.5, 9.2)),))
    im.alpha_composite(band.filter(ImageFilter.GaussianBlur(40)))
    text_center(im, "A dozen tools. A dozen file formats.", 470, font(74, bold=True), TEXT, fade(t, 0.6, 3.8))
    text_center(im, "Weeks of scripts before the first answer.", 470, font(74, bold=True), AMBER, fade(t, 4.0, 6.6))
    text_center(im, "And web servers ask you to upload your data.", 470, font(70, bold=True), TEXT, fade(t, 6.8, 9.4))
    emoji(im, "☁", (W / 2 - 45, 575), 90, fade(t, 7.0, 9.4) * 0.9)
    return im


def s_reveal(t):
    im = background(t + 30).convert("RGBA")
    k = ease(t / 1.2)
    size = int(170 + 120 * k)
    halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(halo).ellipse((W / 2 - size, 360 - size, W / 2 + size, 360 + size), fill=MINT + (int(80 * k),))
    im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(70)))
    s = int(size * 1.15)
    pic = img(str(ICON)).resize((s, s), Image.LANCZOS)
    m = Image.new("L", (s, s), 0)
    ImageDraw.Draw(m).ellipse((0, 0, s, s), fill=int(255 * k))
    pic.putalpha(m)
    im.alpha_composite(pic, (int(W / 2 - s / 2), int(360 - s / 2)))
    text_center(im, "LYNXgwas", 560, font(140, bold=True), TEXT, ease((t - 0.8) / 0.7))
    gradient_text(im, "From summary statistics to the causal variant.", 735, font(56, bold=True), ease((t - 1.7) / 0.7))
    text_center(im, "Every locus, every method, one app on your computer · free and open source", 830, font(36), MUTED,
                ease((t - 2.6) / 0.7))
    return im


def s_steps(t):
    """How easy it is: point it at a file, map the columns, explore."""
    im = background(t + 40).convert("RGBA")
    if t < 3.0:
        section(im, t, "EASY TO USE", "Three steps", "no scripts, no command line", dur=3.0)
        return im
    t -= 3.0
    steps = [("wizard_1.png", (560, 20, 1360, 480), "1", "Name your project", "a short wizard, in the browser"),
             ("wizard_3.png", (560, 20, 1360, 760), "1", "Point it at your summary statistics", "it reads the header and maps the columns"),
             ("home.png", (0, 150, 1920, 700), "2", "It finds and processes the loci", "PLINK clumping, rsIDs, LD, annotation: automatically"),
             ("viewer_tcf7l2.png", (320, 40, 1920, 1080), "3", "Explore every locus", "plots, fine-mapping and genes, one click away")]
    seg = 3.6
    for i, (shot, crop0, num, title, sub) in enumerate(steps):
        a = fade(t, i * seg, (i + 1) * seg + 0.3, 0.45)
        if a <= 0:
            continue
        k = (t - i * seg) / seg
        cw, ch = crop0[2] - crop0[0], crop0[3] - crop0[1]
        crop = lerp_box(crop0, (crop0[0] + cw * 0.02, crop0[1] + ch * 0.02, crop0[2] - cw * 0.02, crop0[3] - ch * 0.02), k)
        frame(im, shot, crop, (620, 70, 1220, 760), a, glow=TEAL)
        halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
        ImageDraw.Draw(halo).ellipse((150, 230, 390, 470), fill=MINT + (int(60 * a),))
        im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(40)))
        text_center(im, num, 205, font(200, bold=True), MINT, a, x=270)
        for j, ln in enumerate(wrap(title, font(50, bold=True), 470)):
            text_at(im, ln, (95, 520 + j * 64), font(50, bold=True), TEXT, a)
        text_at(im, sub, (97, 520 + len(wrap(title, font(50, bold=True), 470)) * 64 + 12), font(30), MUTED, a)
    caption(im, "Your GWAS summary statistics in, every analysis out", "", fade(t, 0, 4 * seg, 0.45), y=900)
    return im


# Feature chapters: (kicker, title, subtitle, accent, items). Each item: (shot, crop or None, title, bullets).
V = (320, 40, 1920, 1080)          # the viewer's main panel
CHAPTERS = [
    ("SEE EVERY LOCUS", "The locus viewer", "everything about a locus on one screen", MINT, [
        ("viewer_pcsk9.png", V, "The genome, then the locus", [
            "genome-wide skyline with every locus marked", "regional plot coloured by LD with the lead SNP",
            "genes, and SNP-heritability in the header"]),
        ("viewer_pcsk9_ld.png", (330, 260, 1910, 1080), "Genes, regulation and LD, together", [
            "histone peaks of the disease's own tissue", "pairwise LD triangle with LD blocks",
            "they zoom together; export as PDF or Excel"]),
        ("viewer_tcf7l2_magma.png", V, "Genes, and the likely causal one", [
            "MAGMA ranks the genes of every locus", "nearest genes and the likely causal gene",
            "split, merge or resize loci live, with undo"]),
    ]),
    ("FINE-MAP IT", "Fine-mapping", "three methods at every locus, and a fourth across datasets", SKYBLUE, [
        ("viewer_sort1_susie.png", V, "SuSiE, ABF and COJO", [
            "95% credible sets drawn on the plot", "run on one locus or on all loci at once",
            "an LD-consistency check flags doubtful SNPs"]),
        ("sum_precision.png", None, "How precise are the sets?", [
            "credible-set sizes across the project", "SuSiE and ABF agree on the lead at 85% of T2D loci",
            "where they differ, a second signal is likely"]),
        ("viewer_sort1_joint.png", V, "Joint fine-mapping across datasets", [
            "mvSuSiE over the datasets of a disease", "SORT1: one two-SNP set, PIP 0.94",
            "sample overlap estimated from null SNPs"]),
    ]),
    ("FIND THE CAUSAL SNP", "Causal-SNP ranking", "one score per SNP from all the evidence", AMBER, [
        ("viewer_pcsk9_rank.png", (336, 325, 1290, 600), "Every SNP, ranked, with its reasons", [
            "SuSiE, ABF, COJO, joint PIP, CADD, LD score, histone peaks, MAGMA",
            "plus how the SNP replicates in the other datasets", "PCSK9 R46L (rs11591147): score 0.993"]),
        ("sum_rankconf.png", None, "How confident, locus by locus", [
            "the top SNP's probability of being causal", "a conditional-logit model trained on FinnGen",
            "held-out chromosomes and diseases for testing"]),
        ("sum_topsnps.png", None, "The strongest candidates", [
            "the most confident causal SNPs of the project", "MTNR1B, ADCY5, CMIP in type 2 diabetes",
            "click a row to open the locus"]),
    ]),
    ("HERITABILITY & QC", "How much do the SNPs explain?", "built-in LD-score regression and local heritability", VIOLET, [
        ("sum_heritability.png", None, "SNP-heritability, genome-wide and per locus", [
            "LDSC on the observed and liability scales", "local heritability (HESS) at every locus",
            "T2D: 0.204; its 132 loci explain 41%"]),
        ("sum_qc.png", None, "Quality control that catches mistakes", [
            "λGC, the LDSC intercept and the ratio", "the sample size implied by the standard errors",
            "a wrong N is flagged and can be corrected"]),
    ]),
    ("GENES & REGULATION", "From SNPs to genes", "gene tests and regulatory enrichment", ORANGE, [
        ("sum_magma.png", None, "MAGMA at every locus", [
            "which genes of the locus carry the association", "125 of 125 tested loci have a significant gene",
            "the strongest gene of each locus"]),
        ("sum_regulatory.png", None, "Regulatory enrichment", [
            "do significant SNPs fall in active chromatin?", "Fisher's exact test per histone mark",
            "T2D in pancreas: OR 2.86 at H3K27ac"]),
    ]),
    ("ACROSS DATASETS", "Compare many GWAS at once", "the same loci, through every dataset", RED, [
        ("serpent.png", (90, 240, 1830, 1060), "Locus Serpent Plot", [
            "927 loci through 30 datasets", "red, thick serpents: loci shared across studies",
            "a category analysis of what drives sharing"]),
        ("serpent_pcsk9.png", (820, 0, 1920, 760), "Click a serpent", [
            "every dataset's lead SNP, OR and p", "PCSK9: one lead SNP in four CAD datasets",
            "joint fine-mapping for every shared locus"]),
        ("constellation.png", (405, 92, 1265, 952), "Gene Constellation", [
            "712 genes across 30 datasets, on one ring", "inner: replicated; color: risk or protective",
            "150 co-significance links between genes"]),
        ("search.png", (480, 55, 1460, 925), "Search every dataset", [
            "a gene, an rsID or a position", "SH2B3: 58 hits across CAD and RA datasets",
            "open any hit straight in the viewer"]),
    ]),
]
ITEM, TITLE_T = 5.4, 2.8


def chapter_scene(ch):
    kicker, title, sub, col, items = ch

    def scene(t):
        im = background(t + 50).convert("RGBA")
        if t < TITLE_T:
            section(im, t, kicker, title, sub, dur=TITLE_T)
            return im
        t -= TITLE_T
        for i, (shot, crop0, ititle, bullets) in enumerate(items):
            a = fade(t, i * ITEM, (i + 1) * ITEM + 0.3, 0.4)
            if a <= 0:
                continue
            src = img(shot)
            c0 = crop0 or (0, 0, src.width, src.height)
            cw, ch_ = c0[2] - c0[0], c0[3] - c0[1]
            k = (t - i * ITEM) / ITEM
            crop = lerp_box(c0, (c0[0] + cw * 0.03, c0[1] + ch_ * 0.03, c0[2] - cw * 0.03, c0[3] - ch_ * 0.03), k)
            frame(im, shot, crop, (720, 70, 1140, 940), a, glow=col if shot.startswith("sum_") else TEAL)
            text_at(im, kicker, (90, 150), font(30, bold=True), col, a)
            y = 205
            for ln in wrap(ititle, font(54, bold=True), 580):
                text_at(im, ln, (90, y), font(54, bold=True), TEXT, a)
                y += 68
            y += 30
            for j, b in enumerate(bullets):
                ab = a * ease((t - i * ITEM - 0.5 - j * 0.45) / 0.4)
                layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
                ImageDraw.Draw(layer).ellipse((92, y + 16, 106, y + 30), fill=col + (int(255 * ab),))
                im.alpha_composite(layer)
                for ln in wrap(b, font(32), 540):
                    text_at(im, ln, (126, y), font(32), (226, 232, 230), ab)
                    y += 44
                y += 18
        return im
    scene.__name__ = "s_ch_" + kicker.split()[0].lower()
    return scene


def s_agent(t):
    im = background(t + 120).convert("RGBA")
    a = fade(t, 0, 7.0, 0.5)
    text_center(im, "ASK THE AGENT", 200, font(38, bold=True), MINT, a)
    gradient_text(im, "Or just ask.", 260, font(110, bold=True), a)
    text_center(im, "a built-in AI agent drives the same actions, with a local model by default", 420, font(38), MUTED, a)
    prompts = ["“Create a project from my T2D summary statistics”", "“Run SuSiE on every locus”",
               "“Which loci are still missing MAGMA?”"]
    for i, p in enumerate(prompts):
        ap = a * ease((t - 1.2 - i * 1.0) / 0.5)
        pill(im, p, (W / 2 - (ImageDraw.Draw(im).textlength(p, font=font(38)) + 44) / 2, 540 + i * 100), font(38),
             TEXT, (22, 52, 58), ap)
    text_center(im, "every step it takes is shown; it cannot delete anything", 900, font(32), MUTED, a * ease((t - 4) / 0.5))
    return im


def s_robust(t):
    """Robust and fast: what was measured."""
    im = background(t + 130).convert("RGBA")
    if t < 2.8:
        section(im, t, "ROBUST AND FAST", "Measured, not promised", "on the article's data, on one desktop computer", dur=2.8)
        return im
    t -= 2.8
    cards = [("4 s", "to read a 7.66-million-SNP GWAS", "and 37 minutes to process all 177 of its loci", MINT),
             ("0.1 SE", "built-in LDSC vs the reference ldsc", "h² 0.7949 vs 0.7967 on identical input", SKYBLUE),
             ("9 of 30", "datasets had a wrong sample size", "caught from their own standard errors, and corrected", AMBER),
             ("2,387", "loci across 30 datasets", "processed and analysed end to end, on one machine", VIOLET)]
    seg = 3.3
    for i, (num, line, sub, col) in enumerate(cards):
        a = fade(t, i * seg, (i + 1) * seg, 0.4)
        if a <= 0:
            continue
        halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
        ImageDraw.Draw(halo).ellipse((W / 2 - 380, 230, W / 2 + 380, 560), fill=col + (int(55 * a),))
        im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(80)))
        text_center(im, num, 250, font(190, bold=True), col, a)
        text_center(im, line, 560, font(58, bold=True), TEXT, a)
        text_center(im, sub, 650, font(34), MUTED, a)
    return im


def bars(im, x, y, label, ours, published, a, col):
    """One pair of horizontal bars: LYNXgwas vs the published estimate."""
    text_at(im, label, (x, y), font(36, bold=True), TEXT, a)
    scale = 1500
    for j, (name, v, c) in enumerate((("LYNXgwas", ours, col), ("published", published, (150, 160, 165)))):
        yy = y + 60 + j * 56
        w = scale * v * ease(a)
        layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
        ImageDraw.Draw(layer).rounded_rectangle((x + 190, yy, x + 190 + w, yy + 40), 8, fill=c + (int(230 * a),))
        im.alpha_composite(layer)
        text_at(im, name, (x, yy + 2), font(30), MUTED, a)
        text_at(im, f"{v:.3f}".rstrip("0").rstrip(".") if v < 1 else str(v), (x + 205 + w, yy + 2), font(30, bold=True), c, a)


def s_published(t):
    im = background(t + 145).convert("RGBA")
    if t < 2.8:
        section(im, t, "DOES IT AGREE?", "Matches the published", "heritability, known causal variants, known loci", dur=2.8)
        return im
    t -= 2.8
    a0 = fade(t, 0, 6.2, 0.45)
    if a0 > 0:
        text_center(im, "SNP-heritability (liability scale)", 130, font(56, bold=True), TEXT, a0)
        bars(im, 120, 290, "Schizophrenia (PGC3)", 0.22, 0.24, a0 * ease((t - 0.4) / 1.0), MINT)
        bars(im, 120, 560, "Type 2 diabetes (Xue 2018)", 0.204, 0.196, a0 * ease((t - 1.2) / 1.0), AMBER)
        text_center(im, "consistent within every disease, and in the published range", 880, font(34), MUTED, a0)
    a1 = fade(t, 6.2, 12.6, 0.45)
    if a1 > 0:
        text_center(im, "Known causal variants, recovered", 130, font(56, bold=True), TEXT, a1)
        rows = [("PCSK9", "R46L · rs11591147", "single-SNP set, PIP 1.00, all 5 CAD datasets", RED),
                ("SORT1", "rs12740374", "two-SNP joint set, PIP 0.94", ORANGE),
                ("SH2B3", "R262W · rs3184504", "lead of a joint set across CAD and RA", AMBER)]
        for i, (gene, snp, what, col) in enumerate(rows):
            ar = a1 * ease((t - 6.6 - i * 0.7) / 0.5)
            y = 300 + i * 170
            layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
            ImageDraw.Draw(layer).rounded_rectangle((260, y, 1660, y + 130), 22, fill=(255, 255, 255, int(14 * ar)),
                                                    outline=col + (int(140 * ar),), width=2)
            im.alpha_composite(layer)
            text_at(im, gene, (300, y + 30), font(60, bold=True, mono=False), col, ar)
            text_at(im, snp, (640, y + 22), font(38, bold=True, mono=True), TEXT, ar)
            text_at(im, what, (640, y + 76), font(30), MUTED, ar)
        text_center(im, "and the classic loci: TCF7L2 · BIN1 · CTLA4 · IL23R · CACNA1C · DRD2", 860, font(34), MUTED,
                    a1 * ease((t - 9) / 0.6))
    return im


def s_found(t):
    im = background(t + 160).convert("RGBA")
    if t < 2.8:
        section(im, t, "WHAT IT FOUND", "Beyond one dataset", "what only a multi-dataset tool can show", dur=2.8)
        return im
    t -= 2.8
    cards = [("60.8%", "of loci: the model's top SNP is the causal one", "vs 50.3% for SuSiE alone (FinnGen in-sample-LD truth)", AMBER),
             ("47%", "of joint credible sets are smaller", "than any single dataset's set: replication narrows the candidates", MINT),
             ("71%", "of which loci are shared is explained by disease", "PERMANOVA across the 30 datasets, p = 0.001", SKYBLUE),
             ("OR 4.4", "enrichment in the disease's own tissue", "rheumatoid arthritis loci in immune-cell H3K4me3 peaks", RED)]
    seg = 3.3
    for i, (num, line, sub, col) in enumerate(cards):
        a = fade(t, i * seg, (i + 1) * seg, 0.4)
        if a <= 0:
            continue
        halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
        ImageDraw.Draw(halo).ellipse((W / 2 - 380, 230, W / 2 + 380, 560), fill=col + (int(55 * a),))
        im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(80)))
        text_center(im, num, 250, font(190, bold=True), col, a)
        text_center(im, line, 560, font(58, bold=True), TEXT, a)
        text_center(im, sub, 650, font(34), MUTED, a)
    return im


DISEASES = [("scz-demo", "Schizophrenia"), ("t2d-demo", "Type 2 diabetes"), ("cad-demo", "Coronary artery disease"),
            ("alz-demo", "Alzheimer's disease"), ("ra-demo", "Rheumatoid arthritis"), ("ibd-demo", "Inflammatory bowel disease")]


def s_corpus(t):
    """Everything above was shown on 30 public GWAS: they come with LYNXgwas as its example data."""
    im = background(t + 175).convert("RGBA")
    a = fade(t, 0, 9.6, 0.5)
    text_center(im, "Everything shown here: 30 public GWAS, 6 diseases", 70, font(60, bold=True), TEXT, a)
    pw, ph, gx, gy = 540, 150, 70, 110
    x0 = (W - (3 * pw + 2 * gx)) / 2
    for i, (pid, name) in enumerate(DISEASES):
        k = ease((t - 0.3 - i * 0.3) / 1.4)
        x = x0 + (i % 3) * (pw + gx)
        y = 240 + (i // 3) * (ph + gy)
        manhattan(im, pid, (x, y, pw, ph), progress=k, alpha=a * ease((t - i * 0.3) / 0.4), small=True, label=name, glow=False)
        text_at(im, f"{SKY[pid]['loci']} loci", (x + pw - 120, y - 40), font(28, bold=True), AMBER, a * ease((t - 1.2 - i * 0.3) / 0.4))
    text_center(im, "2,387 loci · every analysis already run · one click to add them to your LYNXgwas", 800, font(38, bold=True),
                MINT, a * ease((t - 3) / 0.6))
    text_center(im, "the example data downloads from Zenodo; your own projects move as one .zip with Export and Import",
                880, font(32), MUTED, a * ease((t - 4) / 0.6))
    return im


def s_online(t):
    im = background(t + 185).convert("RGBA")
    if t < 3.0:
        section(im, t, "NO INSTALL NEEDED", "Use it online", "LYNXgwas Server: the same tool in your browser", dur=3.0)
        return im
    t -= 3.0
    shots = [("server_home.png", (0, 0, 1920, 1080), (180, 40, 1560, 790),
              "Open the public datasets", "every result and plot, read-only, for anyone"),
             ("server_viewer.png", (320, 40, 1920, 1080), (160, 40, 1600, 790),
              "The same viewer in the browser", "CTLA4 in rheumatoid arthritis"),
             ("server_summary.png", (270, 85, 1650, 860), (180, 40, 1560, 790),
              "Create an account to analyse your own GWAS", "private projects kept for 15 days; or host it for your lab")]
    seg = 3.8
    for i, (shot, crop0, box, title, sub) in enumerate(shots):
        a = fade(t, i * seg, (i + 1) * seg + 0.3, 0.45)
        if a <= 0:
            continue
        k = (t - i * seg) / seg
        crop = lerp_box(crop0, (crop0[0] + 60, crop0[1] + 30, crop0[2] - 60, crop0[3] - 40), k)
        frame(im, shot, crop, box, a, glow=TEAL)
        caption(im, title, sub, a, y=880)
    return im


def s_close(t):
    im = background(t + 210).convert("RGBA")
    a = ease(t / 0.8)
    halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(halo).ellipse((W / 2 - 200, 110, W / 2 + 200, 510), fill=MINT + (int(70 * a),))
    im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(70)))
    s = 260
    pic = img(str(ICON)).resize((s, s), Image.LANCZOS)
    m = Image.new("L", (s, s), 0)
    ImageDraw.Draw(m).ellipse((0, 0, s, s), fill=int(255 * a))
    pic.putalpha(m)
    im.alpha_composite(pic, (W // 2 - s // 2, 180))
    gradient_text(im, "Your GWAS, every locus, answered.", 480, font(92, bold=True), ease((t - 0.6) / 0.7))
    text_center(im, "Free and open source. Runs on your computer; your data stays with you.", 620, font(38), MUTED,
                ease((t - 1.3) / 0.7))
    pa = ease((t - 2.0) / 0.7)
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(layer).rounded_rectangle((W / 2 - 330, 715, W / 2 + 330, 805), 20, fill=(0, 0, 0, int(140 * pa)),
                                            outline=(255, 255, 255, int(50 * pa)))
    im.alpha_composite(layer)
    text_center(im, "pip install lynxgwas", 733, font(46, mono=True), MINT, pa)
    text_center(im, "github.com/AlsammanAlsamman/LYNXgwas", 870, font(44, bold=True), TEXT, ease((t - 2.6) / 0.7))
    text_center(im, "online: LYNXgwas Server · example data on Zenodo", 940, font(32), MUTED, ease((t - 3.0) / 0.7))
    return im


SCENES = ([(s_hook, 13.0), (s_problem, 9.8), (s_reveal, 6.6), (s_steps, 3.0 + 4 * 3.6 + 0.4)]
          + [(chapter_scene(ch), TITLE_T + len(ch[4]) * ITEM + 0.3) for ch in CHAPTERS]
          + [(s_agent, 7.0), (s_robust, 2.8 + 4 * 3.3), (s_published, 2.8 + 12.6), (s_found, 2.8 + 4 * 3.3),
             (s_corpus, 9.6), (s_online, 3.0 + 3 * 3.8 + 0.3), (s_close, 10.0)])


def scene_starts():
    out, s = [], 0.0
    for fn, d in SCENES:
        out.append((fn.__name__, s))
        s += d
    return out


def frame_at(t: float) -> Image.Image:
    start = 0.0
    for scene, dur in SCENES:
        if t < start + dur:
            im = scene(t - start)
            break
        start += dur
    else:
        im = s_close(SCENES[-1][1])
    total = sum(d for _, d in SCENES)
    if t > total - 1.5:
        im.alpha_composite(Image.new("RGBA", im.size, (0, 0, 0, int(255 * ease((t - total + 1.5) / 1.5)))))
    if t < 0.8:
        im.alpha_composite(Image.new("RGBA", im.size, (0, 0, 0, int(255 * (1 - ease(t / 0.8))))))
    return im.convert("RGB")


# --------------------------------------------------------------------------- soundtrack
def soundtrack(path: Path, seconds: float, sr: int = 44100) -> None:
    """Original score: a heartbeat under the opening plot, a dark drone, slow minor pads (Dm-Bb-F-C), a riser
    into the reveal with an impact, a pulsing arpeggio through the tour, and a swell at the numbers."""
    n = int(seconds * sr)
    t = np.arange(n) / sr
    out = np.zeros(n)
    starts = dict(scene_starts())
    reveal, numbers, close = starts["s_reveal"], starts["s_published"], starts["s_close"]

    def note(m):
        return 440.0 * 2 ** ((m - 69) / 12)

    def env(start, dur, attack=1.5, release=2.0):
        e = np.clip((t - start) / attack, 0, 1) * np.clip((start + dur - t) / release, 0, 1)
        return e * ((t >= start) & (t <= start + dur))

    # heartbeat (lub-dub) under the hook, speeding up into the problem
    beat, k = 0.6, 0
    while beat < reveal - 0.5:
        for off, amp in ((0.0, 1.0), (0.22, 0.6)):
            s0 = beat + off
            idx = (t >= s0) & (t < s0 + 0.35)
            tt = t[idx] - s0
            out[idx] += 0.32 * amp * np.sin(2 * np.pi * (55 + 30 * np.exp(-tt * 18)) * tt) * np.exp(-tt * 11)
        k += 1
        beat += max(0.55, 1.05 - 0.018 * k)
    out += 0.09 * np.sin(2 * np.pi * note(26) * t) * np.clip(t / 5, 0, 1)                      # D1 drone
    out += 0.04 * np.sin(2 * np.pi * note(38) * 1.002 * t) * np.clip((t - 4) / 6, 0, 1)
    chords = [(62, 65, 69), (58, 62, 65), (53, 57, 60), (60, 64, 67)]                         # Dm Bb F C
    bar = 4.0
    for b in range(int(seconds / bar) + 1):
        s0 = b * bar
        if s0 < 6:
            continue
        for m in chords[b % 4]:
            f = note(m - 12)
            e = env(s0, bar + 1.5, attack=1.4, release=1.8)
            out += 0.032 * e * (np.sin(2 * np.pi * f * t) + 0.5 * np.sin(2 * np.pi * f * 1.004 * t)
                                + 0.22 * np.sin(2 * np.pi * f * 2 * t))
    # riser and impact into the reveal
    r0 = reveal - 5
    rs = (t > r0) & (t < reveal)
    out[rs] += 0.06 * np.sin(2 * np.pi * (180 + (t[rs] - r0) ** 2 * 70) * t[rs]) * ((t[rs] - r0) / 5) ** 2
    noise = np.random.default_rng(1).standard_normal(n)
    out[rs] += 0.03 * noise[rs] * ((t[rs] - r0) / 5) ** 3
    hit = (t >= reveal) & (t < reveal + 5)
    out[hit] += 0.30 * np.sin(2 * np.pi * note(38) * t[hit]) * np.exp(-(t[hit] - reveal) * 1.2)
    out[hit] += 0.18 * np.sin(2 * np.pi * note(26) * t[hit]) * np.exp(-(t[hit] - reveal) * 0.8)
    # pulsing arpeggio from the corpus scene to the close
    arp0, step = starts["s_steps"], 0.25
    for j in range(int((close - arp0) / step)):
        s0 = arp0 + j * step
        b = int(s0 / bar) % 4
        m = chords[b][j % 3] + (12 if j % 8 >= 6 else 0)
        idx = (t >= s0) & (t < s0 + 1.0)
        tt = t[idx] - s0
        out[idx] += 0.024 * np.sin(2 * np.pi * note(m) * tt) * np.exp(-tt * 5)
        if j % 4 == 0:                                                                          # soft pulse
            idx2 = (t >= s0) & (t < s0 + 0.25)
            t2 = t[idx2] - s0
            out[idx2] += 0.10 * np.sin(2 * np.pi * 60 * t2) * np.exp(-t2 * 16)
    sw = (t > numbers - 2) & (t < numbers + 14)
    out[sw] += 0.05 * np.sin(2 * np.pi * note(50) * t[sw]) * np.sin(np.pi * (t[sw] - numbers + 2) / 16) ** 2
    end = (t >= close) & (t < seconds)
    for m in (50, 57, 62, 65):
        out[end] += 0.035 * np.sin(2 * np.pi * note(m) * t[end]) * np.clip((t[end] - close) / 1.5, 0, 1)
    wet = np.zeros(n)
    for delay, g in ((0.11, 0.35), (0.23, 0.25), (0.41, 0.16), (0.67, 0.08)):
        d = int(delay * sr)
        wet[d:] += g * out[:-d]
    out = out + wet
    out *= np.clip((seconds - t) / 3.0, 0, 1)
    out /= max(1e-9, np.abs(out).max()) / 0.82
    stereo = np.stack([out, np.roll(out, int(0.013 * sr))], axis=1)
    with wave.open(str(path), "wb") as wf:
        wf.setnchannels(2)
        wf.setsampwidth(2)
        wf.setframerate(sr)
        wf.writeframes((stereo * 32767).astype(np.int16).tobytes())


# --------------------------------------------------------------------------- render
def main() -> None:
    total = sum(d for _, d in SCENES)
    if "--preview" in sys.argv:
        secs = [float(x) for x in sys.argv[2:]] or [3, 12, 20, 30, 40, 52, 66, 78, 92, 104, 118, 130, 140, 150,
                                                   162, 172, 184, 196]
        for sec in secs:
            frame_at(sec).save(HERE / f"preview_{int(sec):03d}.png")
        print("preview frames written; scenes:", [(n, round(s, 1)) for n, s in scene_starts()], f"total {total:.1f}s")
        return
    import imageio_ffmpeg
    audio = HERE / "soundtrack.wav"
    soundtrack(audio, total)
    video = HERE / "lynxgwas-intro.mp4"
    proc = subprocess.Popen([
        imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error",
        "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-",
        "-i", str(audio),
        "-c:v", "libx264", "-preset", "medium", "-crf", "20", "-pix_fmt", "yuv420p",
        "-c:a", "aac", "-b:a", "192k", "-shortest", "-movflags", "+faststart", str(video)],
        stdin=subprocess.PIPE)
    frames = int(total * FPS)
    for i in range(frames):
        proc.stdin.write(frame_at(i / FPS).tobytes())
        if i % (FPS * 10) == 0:
            print(f"  {i / FPS:5.0f}s / {total:.0f}s", flush=True)
    proc.stdin.close()
    proc.wait()
    print(f"Wrote {video} ({video.stat().st_size / 1e6:.1f} MB, {total:.0f} s)")


if __name__ == "__main__":
    main()
