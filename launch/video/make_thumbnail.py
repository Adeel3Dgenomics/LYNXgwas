"""YouTube thumbnail (1280x720) and the README's clickable video image (with a play button).

    python launch/video/make_thumbnail.py
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).parent))
from make_video import (AMBER, ICON, MINT, MUTED, TEXT, background, font, frame, img,  # noqa: E402
                        manhattan)

ROOT = Path(__file__).resolve().parents[2]


def thumbnail(play_button: bool) -> Image.Image:
    im = background(3).convert("RGBA")
    manhattan(im, "scz-demo", (60, 790, 1800, 230), alpha=0.55, glow=True)
    frame(im, "viewer_pcsk9_rank.png", (336, 225, 1500, 620), (1010, 230, 860, 520), 1.0, glow=(120, 80, 20))
    s = 190
    pic = img(str(ICON)).resize((s, s), Image.LANCZOS)
    m = Image.new("L", (s, s), 0)
    ImageDraw.Draw(m).ellipse((0, 0, s, s), fill=255)
    pic.putalpha(m)
    im.alpha_composite(pic, (100, 90))
    d = ImageDraw.Draw(im)
    d.text((320, 110), "LYNXgwas", font=font(120, bold=True), fill=TEXT)
    y = 360
    for line, col in (("Find the causal", TEXT), ("variant.", MINT)):
        d.text((100, y), line, font=font(112, bold=True), fill=col)
        y += 128
    d.text((104, y + 20), "30 GWAS · 6 diseases · one app", font=font(54, bold=True), fill=AMBER)
    d.text((104, y + 96), "free, open source, runs on your computer", font=font(40), fill=MUTED)
    if play_button:
        ring = Image.new("RGBA", im.size, (0, 0, 0, 0))
        rd = ImageDraw.Draw(ring)
        cx, cy, r = 1440, 490, 115
        rd.ellipse((cx - r, cy - r, cx + r, cy + r), fill=(0, 0, 0, 175))
        rd.polygon([(cx - 36, cy - 60), (cx - 36, cy + 60), (cx + 64, cy)], fill=(255, 255, 255, 240))
        im.alpha_composite(ring)
    return im.convert("RGB").resize((1280, 720), Image.LANCZOS)


if __name__ == "__main__":
    thumbnail(False).save(Path(__file__).parent / "youtube-thumbnail.png")
    thumbnail(True).save(ROOT / "docs" / "images" / "video-thumbnail.png")
    print("thumbnails written")
