"""Render the LYNXgwas intro video (1080p, about 3 minutes) with an original synthesized soundtrack.

    pip install pillow numpy imageio-ffmpeg
    python launch/video/capture.py && python launch/video/capture_extra.py   # screenshots (app running)
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
    """One GWAS draws itself: 7.66 million SNPs, and then the question."""
    im = background(t).convert("RGBA")
    box = (110, 260, W - 220, 520)
    manhattan(im, "scz-demo", box, progress=ease(t / 6.5), alpha=ease(t / 0.8) * (1 - 0.55 * ease((t - 9) / 1.2)))
    n = int(SKY["scz-demo"]["snps"] * ease(t / 6.5))
    text_center(im, f"{n:,} SNPs", 120, font(80, bold=True, mono=True), TEXT, fade(t, 0.3, 7.6))
    text_center(im, "one schizophrenia GWAS (PGC3)", 225, font(34), MUTED, fade(t, 0.8, 7.6))
    text_center(im, "177 genome-wide significant loci.", 120, font(76, bold=True), AMBER, fade(t, 7.8, 10.8))
    q = [("Which variant is causal?", 9.2), ("Which gene does it act through?", 10.4), ("Does it replicate in other studies?", 11.6)]
    for i, (line, at) in enumerate(q):
        text_center(im, line, 830 + i * 64 - 70, font(50, bold=True), TEXT if i < 2 else MINT, fade(t, at, 15.0))
    return im


TOOLS = ["PLINK", "SuSiE", "FINEMAP / ABF", "GCTA-COJO", "MAGMA", "LDSC", "HESS", "mvSuSiE", "coloc", "SuSiEx",
         "GWAMA", "liftOver", "LocusZoom", "R scripts", "bash", "Excel", "reference panels", "GENCODE",
         "Roadmap peaks", "HapMap3", "CADD", "1000 Genomes", "python", "awk"]


def s_problem(t):
    """A dozen tools and formats swirling; then: and the web servers want your data."""
    im = background(t + 15).convert("RGBA")
    rnd = random.Random(3)
    collapse = ease((t - 10.2) / 1.6)
    for i, name in enumerate(TOOLS):
        x0 = rnd.uniform(140, W - 360)
        y0 = rnd.choice((rnd.uniform(90, 330), rnd.uniform(700, H - 120)))   # keep the headline band clear
        vx, vy = rnd.uniform(-22, 22), rnd.uniform(-6, 6)
        x = x0 + vx * t + math.sin(t * 1.3 + i) * 18
        y = y0 + vy * t + math.cos(t * 1.1 + i) * 14
        x = x + (W / 2 - 80 - x) * collapse
        y = y + (H / 2 - y) * collapse
        a = ease((t - i * 0.12) / 0.5) * (1 - collapse) * 0.9
        size = rnd.choice((34, 40, 46, 52))
        col = rnd.choice((TEXT, MUTED, SKYBLUE, AMBER, MINT))
        text_at(im, name, (x, y), font(size, bold=size > 40, mono=rnd.random() < 0.3), col, a)
    band = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(band).rectangle((0, 400, W, 690), fill=BG + (int(240 * fade(t, 0.5, 10.4)),))
    im.alpha_composite(band.filter(ImageFilter.GaussianBlur(40)))
    text_center(im, "A dozen tools. A dozen file formats.", 470, font(74, bold=True), TEXT, fade(t, 0.6, 4.6))
    text_center(im, "Weeks of scripts before the first answer.", 470, font(74, bold=True), AMBER, fade(t, 4.8, 7.6))
    text_center(im, "And web servers ask you to upload your data.", 470, font(70, bold=True), TEXT, fade(t, 7.8, 10.6))
    emoji(im, "☁", (W / 2 - 55, 575), 90, fade(t, 8.2, 10.6) * 0.9)
    return im


def s_reveal(t):
    im = background(t + 30).convert("RGBA")
    k = ease(t / 1.2)
    size = int(170 + 120 * k)
    halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(halo).ellipse((W / 2 - size, 360 - size, W / 2 + size, 360 + size), fill=MINT + (int(80 * k),))
    im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(70)))
    ic = img(str(ICON))
    s = int(size * 1.15)
    pic = ic.resize((s, s), Image.LANCZOS)
    m = Image.new("L", (s, s), 0)
    ImageDraw.Draw(m).ellipse((0, 0, s, s), fill=int(255 * k))
    pic.putalpha(m)
    im.alpha_composite(pic, (int(W / 2 - s / 2), int(360 - s / 2)))
    text_center(im, "LYNXgwas", 560, font(140, bold=True), TEXT, ease((t - 0.8) / 0.7))
    gradient_text(im, "Every locus. Every method. One app on your computer.", 735, font(54, bold=True),
                  ease((t - 1.7) / 0.7))
    text_center(im, "Locus analYsis and geNomic eXplorer · free and open source", 830, font(36), MUTED,
                ease((t - 2.6) / 0.7))
    return im


DISEASES = [("scz-demo", "Schizophrenia"), ("t2d-demo", "Type 2 diabetes"), ("cad-demo", "Coronary artery disease"),
            ("alz-demo", "Alzheimer's disease"), ("ra-demo", "Rheumatoid arthritis"), ("ibd-demo", "Inflammatory bowel disease")]


def s_corpus(t):
    """The test: 30 public GWAS, 6 diseases. Six Manhattan plots fill the screen, then the counters."""
    im = background(t + 40).convert("RGBA")
    if t < 3.4:
        section(im, t, "THE TEST", "30 public GWAS", "six diseases, five datasets each, every locus processed")
        return im
    t -= 3.4
    a = ease(t / 0.5) * ease((11.6 - t) / 0.5)
    pw, ph, gx, gy = 540, 160, 70, 120
    x0 = (W - (3 * pw + 2 * gx)) / 2
    for i, (pid, name) in enumerate(DISEASES):
        k = ease((t - i * 0.35) / 1.6)
        x = x0 + (i % 3) * (pw + gx)
        y = 150 + (i // 3) * (ph + gy)
        manhattan(im, pid, (x, y, pw, ph), progress=k, alpha=a * ease((t - i * 0.35) / 0.4), small=True,
                  label=name, glow=False)
        text_at(im, f"{SKY[pid]['loci']} loci", (x + pw - 120, y - 40), font(28, bold=True), AMBER,
                a * ease((t - 1.4 - i * 0.35) / 0.4))
    stats = [("30", "datasets"), ("6", "diseases"), ("2,387", "loci"), ("927", "loci aligned across them")]
    for i, (num, lab) in enumerate(stats):
        k = ease((t - 4.0 - i * 0.7) / 0.6) * a
        cx = W * (i + 0.5) / 4
        text_center(im, num, 770, font(96, bold=True), MINT if i < 3 else AMBER, k, x=cx)
        text_center(im, lab, 890, font(32), MUTED, k, x=cx)
    text_center(im, "all on GRCh37, all processed on one computer", 975, font(30), MUTED, ease((t - 7) / 0.6) * a)
    return im


PILLARS = [
    ("🎯", "Fine-mapping", "SuSiE, Wakefield ABF and GCTA-COJO at every locus", MINT),
    ("🔗", "Joint fine-mapping", "mvSuSiE across the datasets of a disease", SKYBLUE),
    ("🏆", "Causal-SNP ranking", "one model combining all the evidence", AMBER),
    ("🧬", "Genes & heritability", "MAGMA per locus, LDSC and local HESS", VIOLET),
    ("🗺️", "Across datasets", "Locus Serpent Plot and Gene Constellation", ORANGE),
    ("🔒", "Private by design", "runs on your computer; nothing is uploaded", RED),
]


def s_pillars(t):
    im = background(t + 60).convert("RGBA")
    a = ease(t / 0.6) * ease((15.4 - t) / 0.6)
    text_center(im, "From summary statistics to the causal variant", 85, font(70, bold=True), TEXT, a)
    text_center(im, "every method, on every locus, in one place", 185, font(36), MUTED, a)
    tw, th, gx, gy = 540, 270, 40, 40
    x0 = (W - (3 * tw + 2 * gx)) / 2
    for i, (em, name, line, col) in enumerate(PILLARS):
        ap = ease((t - 1.0 - i * 1.3) / 0.6) * a
        if ap <= 0:
            continue
        x = x0 + (i % 3) * (tw + gx)
        y = 290 + (i // 3) * (th + gy) + (1 - ap) * 40
        card = Image.new("RGBA", im.size, (0, 0, 0, 0))
        cd = ImageDraw.Draw(card)
        cd.rounded_rectangle((x, y, x + tw, y + th), 26, fill=(255, 255, 255, int(16 * ap)),
                             outline=col + (int(110 * ap),), width=2)
        im.alpha_composite(card)
        emoji(im, em, (x + 34, y + 40), 62, ap)
        text_at(im, name, (x + 112, y + 48), font(40, bold=True), col, ap)
        for j, ln in enumerate(wrap(line, font(33), tw - 70)):
            text_at(im, ln, (x + 36, y + 138 + j * 46), font(33), (226, 232, 230), ap)
    return im


def s_home(t):
    im = background(t + 80).convert("RGBA")
    a = fade(t, 0, 8.0, 0.5)
    box = lerp_box((160, 60, 1600, 760), (110, 40, 1700, 800), t / 8)
    frame(im, "home.png", (0, 0, 1920, 1080), box, a, glow=TEAL)
    caption(im, "All 30 datasets, one home page", "each one processed, annotated and up to date", a, y=890)
    return im


def s_pcsk9(t):
    """The showcase: PCSK9 R46L, every method agrees, the ranking finds it, the LD triangle shows why."""
    im = background(t + 90).convert("RGBA")
    if t < 3.4:
        section(im, t, "A LOCUS, CLOSE UP", "PCSK9", "coronary artery disease · chr1:55.5 Mb")
        return im
    t -= 3.4
    seg = [(0, 6.5), (6.5, 12.5), (12.5, 18.6)]
    a0 = fade(t, *seg[0], 0.5)
    if a0 > 0:
        crop = lerp_box((0, 0, 1920, 1080), (330, 200, 1910, 800), (t - 0.8) / 3.5)
        frame(im, "viewer_pcsk9.png", crop, (100, 40, 1720, 780), a0, glow=TEAL)
        caption(im, "Every method agrees on rs11591147", "SuSiE PIP 1.00 · ABF 1.00 · COJO selected · joint mvSuSiE 1.00 in 5 CAD datasets",
                a0, y=880)
    a1 = fade(t, *seg[1], 0.5)
    if a1 > 0:
        crop = lerp_box((336, 230, 1910, 640), (336, 325, 1340, 560), (t - 6.8) / 4.5)
        frame(im, "viewer_pcsk9_rank.png", crop, (100, 80, 1720, 700), a1, glow=(120, 80, 20))
        caption(im, "The causal-SNP ranking puts it first: score 0.993",
                "PCSK9 R46L, a known causal missense variant, found from summary statistics alone", a1, y=880)
    a2 = fade(t, *seg[2], 0.5)
    if a2 > 0:
        crop = lerp_box((330, 560, 1910, 1080), (600, 760, 1700, 1080), (t - 12.8) / 5)
        frame(im, "viewer_pcsk9_ld.png", crop, (100, 60, 1720, 740), a2, glow=TEAL)
        caption(im, "Regulatory tracks, genes and LD, zooming together",
                "histone peaks of the disease tissue, and the LD triangle of the locus", a2, y=880)
    return im


MONTAGE = [("viewer_tcf7l2.png", "TCF7L2 · type 2 diabetes", "the strongest T2D signal; lead SNP rs7903146"),
           ("viewer_il23r.png", "IL23R · Crohn's disease", "p = 6 × 10⁻⁹⁷"),
           ("viewer_drd2.png", "DRD2 · schizophrenia", "the dopamine receptor locus"),
           ("viewer_sort1.png", "SORT1 · coronary artery disease", "rs12740374, the LDL-cholesterol variant")]


def s_montage(t):
    im = background(t + 110).convert("RGBA")
    seg = 3.4
    for i, (shot, title, sub) in enumerate(MONTAGE):
        a = fade(t, i * seg, (i + 1) * seg + 0.3, 0.45)
        if a <= 0:
            continue
        k = (t - i * seg) / seg
        crop = lerp_box((320, 40, 1920, 1080), (330, 200, 1910, 820), k)
        frame(im, shot, crop, (160, 50, 1600, 760), a, glow=TEAL)
        caption(im, title, sub, a, y=880)
    return im


def s_serpent(t):
    im = background(t + 125).convert("RGBA")
    if t < 3.2:
        section(im, t, "ACROSS ALL 30 DATASETS", "Locus Serpent Plot", "927 loci, each one a serpent through every dataset", dur=3.2)
        return im
    t -= 3.2
    a0 = fade(t, 0, 6.8, 0.5)
    if a0 > 0:
        crop = lerp_box((90, 240, 1830, 1060), (380, 280, 1250, 880), (t - 1.0) / 5)
        frame(im, "serpent.png", crop, (90, 40, 1740, 790), a0, glow=(120, 30, 30))
        caption(im, "Shared loci light up as long red serpents",
                "color and thickness = significance in each dataset; drift = where its lead SNP sits", a0, y=880)
    a1 = fade(t, 6.8, 12.2, 0.5)
    if a1 > 0:
        crop = lerp_box((820, 0, 1920, 760), (830, 0, 1920, 560), (t - 7) / 5)
        frame(im, "serpent_pcsk9.png", crop, (240, 40, 1440, 790), a1, glow=(120, 30, 30))
        caption(im, "Click a serpent: every dataset's lead SNP, OR and p",
                "PCSK9: the same lead SNP in four coronary artery disease datasets", a1, y=880)
    return im


def s_constellation(t):
    im = background(t + 140).convert("RGBA")
    a = fade(t, 0, 9.0, 0.6)
    s = img("constellation.png").crop((405, 92, 1265, 952))
    zoom = 1 + 0.10 * ease(t / 9)
    size = int(760 * zoom)
    rot = s.rotate(-6 * t, resample=Image.BICUBIC, fillcolor=(255, 255, 255, 255)).resize((size, size), Image.BILINEAR)
    m = Image.new("L", (size, size), 0)
    ImageDraw.Draw(m).ellipse((8, 8, size - 8, size - 8), fill=int(255 * a))
    rot.putalpha(m)
    halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
    ImageDraw.Draw(halo).ellipse((560 - size / 2 - 40, 450 - size / 2 - 40, 560 + size / 2 + 40, 450 + size / 2 + 40),
                                 fill=VIOLET + (int(70 * a),))
    im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(60)))
    im.alpha_composite(rot, (int(560 - size / 2), int(450 - size / 2)))
    text_at(im, "Gene Constellation", (1050, 290), font(80, bold=True), TEXT, a)
    for j, ln in enumerate(["712 genes from 30 datasets", "inner ring: replicated in more datasets",
                            "color: risk or protective", "150 co-significance links"]):
        text_at(im, ln, (1056, 420 + j * 62), font(38), MUTED if j else AMBER, ease((t - 0.8 - j * 0.6) / 0.6) * a)
    return im


def s_summary(t):
    """A project's summary page, scrolling: every analysis on every locus."""
    im = background(t + 150).convert("RGBA")
    a = fade(t, 0, 12.0, 0.5)
    src = img("summary_t2d.png")
    y0 = (src.height - 1080) * ease((t - 1.0) / 9.5) * 0.72
    frame(im, "summary_t2d.png", (60, y0, 1860, y0 + 1000), (180, 40, 1560, 790), a, glow=TEAL)
    caption(im, "Every analysis, every locus, at a glance",
            "progress, heritability, credible sets, ranking, MAGMA and what still needs to run", a, y=880)
    return im


NUMBERS = [("60.8%", "of loci: the model's top SNP is the causal one", "vs 50.3% for SuSiE alone · FinnGen in-sample-LD truth", AMBER),
           ("47%", "of joint credible sets are smaller", "than any single dataset's, from mvSuSiE across datasets", MINT),
           ("3,226", "credible sets from 2,299 fine-mapped loci", "SuSiE, ABF and COJO agree on the lead SNP at most loci", SKYBLUE),
           ("0 bytes", "uploaded anywhere", "every analysis ran on one desktop computer", RED)]


def s_numbers(t):
    im = background(t + 165).convert("RGBA")
    seg = 3.6
    for i, (num, line, sub, col) in enumerate(NUMBERS):
        a = fade(t, i * seg, (i + 1) * seg, 0.45)
        if a <= 0:
            continue
        k = ease((t - i * seg) / 1.2)
        shown = num
        digits = "".join(c for c in num if c.isdigit() or c == ".")
        if digits and num != "0 bytes":
            val = float(digits.replace(",", "")) * k
            shown = (f"{val:.1f}%" if "." in num else (f"{int(val):,}" + ("%" if "%" in num else "")))
        halo = Image.new("RGBA", im.size, (0, 0, 0, 0))
        ImageDraw.Draw(halo).ellipse((W / 2 - 380, 230, W / 2 + 380, 560), fill=col + (int(55 * a),))
        im.alpha_composite(halo.filter(ImageFilter.GaussianBlur(80)))
        text_center(im, shown, 250, font(200, bold=True), col, a)
        text_center(im, line, 560, font(58, bold=True), TEXT, a)
        text_center(im, sub, 650, font(34), MUTED, a)
    text_center(im, "Results of the LYNXgwas validation on the 30-dataset corpus", 960, font(28), MUTED, fade(t, 0, 4 * seg))
    return im


def s_online(t):
    im = background(t + 180).convert("RGBA")
    if t < 3.4:
        section(im, t, "ALSO ONLINE", "LYNXgwas Server", "the same app for a whole lab, or the world")
        return im
    t -= 3.4
    shots = [("server_home.png", (0, 0, 1920, 1080), (180, 40, 1560, 790),
              "Public datasets, open to everyone", "read-only, with every result and plot"),
             ("server_viewer.png", (320, 40, 1920, 1080), (160, 40, 1600, 790),
              "The same viewer in the browser", "CTLA4 in rheumatoid arthritis, served by LYNXgwas Server"),
             ("server_summary.png", (330, 90, 1590, 830), (180, 40, 1560, 790),
              "Every summary, open in the browser", "visitors browse read-only; an account analyses its own data")]
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


def s_data(t):
    im = background(t + 195).convert("RGBA")
    a0 = fade(t, 0, 6.6, 0.5)
    if a0 > 0:
        crop = lerp_box((0, 150, 1920, 520), (130, 160, 1800, 480), t / 6.6)
        frame(im, "home_example.png", crop, (120, 160, 1680, 520), a0, glow=TEAL)
        caption(im, "Start with the article's data in one click", "30 datasets with every result, downloaded from Zenodo and added for you",
                a0, y=780)
    a1 = fade(t, 6.6, 13.0, 0.5)
    if a1 > 0:
        crop = lerp_box((630, 40, 1290, 670), (640, 60, 1280, 650), (t - 6.6) / 6.4)
        frame(im, "export_dialog.png", crop, (560, 40, 800, 760), a1, glow=TEAL)
        caption(im, "Move projects as one .zip", "Export one dataset or many; Import them on any computer", a1, y=880)
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
    text_center(im, "example data on Zenodo · LYNXgwas Server for your lab", 940, font(32), MUTED, ease((t - 3.0) / 0.7))
    return im


SCENES = [(s_hook, 15.0), (s_problem, 11.4), (s_reveal, 7.0), (s_corpus, 15.4), (s_pillars, 16.0), (s_home, 8.0),
          (s_pcsk9, 22.0), (s_montage, 13.9), (s_serpent, 15.4), (s_constellation, 9.0), (s_summary, 12.0),
          (s_numbers, 14.4), (s_online, 15.2), (s_data, 13.0), (s_close, 10.0)]


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
    reveal, numbers, close = starts["s_reveal"], starts["s_numbers"], starts["s_close"]

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
    arp0, step = starts["s_corpus"], 0.25
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
