"""Feature screenshots for the intro video: project wizard, search, credible sets on the plot, summary sections.

    python launch/video/capture_features.py      (desktop app running on http://localhost:8765)
"""
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

OUT = Path(__file__).parent / "shots"
DESK = "http://127.0.0.1:8765"
SECTIONS = {"sum_heritability": "Heritability", "sum_qc": "Quality control",
            "sum_signals": "Fine-mapping: how many signals per locus", "sum_precision": "Fine-mapping: how precise are the credible sets",
            "sum_rankconf": "Causal-SNP ranking: confidence per locus", "sum_topsnps": "Top causal-SNP candidates",
            "sum_magma": "Genes (MAGMA)", "sum_regulatory": "Regulatory enrichment", "sum_progress": "Analysis progress",
            "sum_matrix": "Which loci have which results", "sum_joint": "Joint fine-mapping across datasets (mvSuSiE)",
            "sum_sce2g": "Enhancer–gene links (scE2G)"}


def save(page, name):
    page.screenshot(path=str(OUT / f"{name}.png"))
    print("  ", name, flush=True)


def main():
    with sync_playwright() as p:
        b = p.chromium.launch(channel="msedge", headless=True)
        page = b.new_page(viewport={"width": 1920, "height": 1080})

        # the project wizard, on an existing project so every step is filled in
        page.goto(f"{DESK}/", wait_until="domcontentloaded")
        time.sleep(5)
        page.evaluate("editConfig('t2d-demo')")
        time.sleep(3)
        page.fill("[data-wiz-key=_name]", "Type 2 diabetes (Xue et al. 2018)")
        page.fill("[data-wiz-key=_desc]", "T2D GWAS, 62,892 cases and 596,424 controls")
        save(page, "wizard_1")
        for i in range(2, 5):
            try:
                page.click("#wiz-next", timeout=5000)
                time.sleep(2.5)
                if i != 2:   # step 2 shows the local file paths (user name, folders): never saved
                    save(page, f"wizard_{i}")
            except Exception as e:
                print("   wizard step", i, e)
                break

        # search every dataset at once
        page.goto(f"{DESK}/", wait_until="domcontentloaded")
        time.sleep(5)
        page.evaluate("openGlobalSearch()")
        time.sleep(2)
        box = page.locator("#gsearch-box")
        box.locator("input[type=text], input:not([type])").first.fill("SH2B3")
        box.get_by_role("button", name="Search").first.click()
        time.sleep(8)
        save(page, "search")

        # credible sets drawn on the locus plot (SORT1: a two-SNP joint set at PIP 0.94)
        page.goto(f"{DESK}/viewer.html?project=cad-demo&locus=3", wait_until="domcontentloaded", timeout=180_000)
        time.sleep(14)
        chips = page.get_by_text("Show on plot", exact=True)
        try:
            chips.nth(0).click()
            time.sleep(4)
            save(page, "viewer_sort1_susie")
            chips.nth(chips.count() - 1).click()
            time.sleep(5)
            save(page, "viewer_sort1_joint")
        except Exception as e:
            print("   show on plot:", e)

        # the summary page, one section at a time
        page.close()
        page = b.new_page(viewport={"width": 1600, "height": 1000}, device_scale_factor=2)
        page.goto(f"{DESK}/summary.html?project=t2d-demo", wait_until="domcontentloaded")
        time.sleep(30)
        for name, title in SECTIONS.items():
            try:
                sec = page.locator("section.card", has=page.locator("h2", has_text=title)).first
                sec.scroll_into_view_if_needed()
                time.sleep(1.5)
                sec.screenshot(path=str(OUT / f"{name}.png"))
                if name == "sum_sce2g":   # a white margin, so the video's slow zoom on its top part clips nothing
                    from PIL import Image
                    im = Image.open(OUT / f"{name}.png")
                    pad = Image.new("RGB", (im.width + 160, im.height + 80), "white")
                    pad.paste(im, (80, 80))
                    pad.save(OUT / f"{name}.png")
                print("  ", name, flush=True)
            except Exception as e:
                print("   section", name, e)
        b.close()


if __name__ == "__main__":
    main()
