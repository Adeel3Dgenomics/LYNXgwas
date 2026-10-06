"""Screenshots for the intro video, taken from a running LYNXgwas with the article's 30 datasets.

    desktop app on http://localhost:8765, server edition on http://localhost:8780
    python launch/video/capture.py      ->  launch/video/shots/*.png

Uses Playwright with the system's Microsoft Edge (no browser download).
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

OUT = Path(__file__).parent / "shots"
DESK = "http://localhost:8765"
SERVER = "http://localhost:8780"
SERPENT = "d5924980-27fb-497a-8443-1ebc69d076e1"
CONSTELLATION = "ef832895-6854-454d-9265-bcd5f7acf9f5"


def shot(page, name, url=None, wait=6.0, full=False, before=None, selector_wait=None):
    if url:
        page.goto(url, wait_until="domcontentloaded", timeout=180_000)
    if selector_wait:
        try:
            page.wait_for_selector(selector_wait, timeout=120_000)
        except Exception as e:
            print("  (selector wait failed:", e, ")")
    time.sleep(wait)
    if before:
        before(page)
        time.sleep(3)
    page.screenshot(path=str(OUT / f"{name}.png"), full_page=full)
    print("  ", name, flush=True)


def main(only=None):
    OUT.mkdir(exist_ok=True)
    with sync_playwright() as p:
        b = p.chromium.launch(channel="msedge", headless=True)
        page = b.new_page(viewport={"width": 1920, "height": 1080}, device_scale_factor=1)

        jobs = [
            ("home", f"{DESK}/", 5, False, None, None),
            ("home_example", f"{DESK}/", 5, False, lambda pg: pg.evaluate("exdataShow()"), None),
            ("export_dialog", f"{DESK}/", 5, False,
             lambda pg: (pg.evaluate("packOpenExport()"), pg.evaluate(
                 "document.querySelectorAll('.pack-pick').forEach((c,i)=>{if(i<9)c.checked=true})")), None),
            ("agent", f"{DESK}/", 5, False, lambda pg: pg.evaluate("agentToggle()"), None),
            ("viewer_pcsk9", f"{DESK}/viewer.html?project=cad-demo&locus=1", 14, False, None, None),
            ("viewer_sort1", f"{DESK}/viewer.html?project=cad-demo&locus=3", 14, False, None, None),
            ("viewer_tcf7l2", f"{DESK}/viewer.html?project=t2d-demo&locus=87", 14, False, None, None),
            ("viewer_il23r", f"{DESK}/viewer.html?project=ibd-demo&locus=2", 14, False, None, None),
            ("viewer_drd2", f"{DESK}/viewer.html?project=scz-demo&locus=118", 14, False, None, None),
            ("viewer_tcf7l2_full", f"{DESK}/viewer.html?project=t2d-demo&locus=87", 14, True, None, None),
            ("summary_t2d", f"{DESK}/summary.html?project=t2d-demo", 25, True, None, None),
            ("summary_cad", f"{DESK}/summary.html?project=cad-demo", 25, True, None, None),
            ("serpent", f"{DESK}/serpent_plot.html?job={SERPENT}", 25, False, None, None),
            ("serpent_full", f"{DESK}/serpent_plot.html?job={SERPENT}", 5, True, None, None),
            ("serpent_pcsk9", f"{DESK}/serpent_plot.html?job={SERPENT}&locus=27", 25, False, None, None),
            ("constellation", f"{DESK}/gene_constellation.html?job={CONSTELLATION}", 30, False, None, None),
            ("constellation_full", f"{DESK}/gene_constellation.html?job={CONSTELLATION}", 5, True, None, None),
            ("server_home", f"{SERVER}/", 8, False, None, None),
            ("server_viewer", f"{SERVER}/viewer.html?project=ra-demo&locus=16", 16, False, None, None),
            ("server_summary", f"{SERVER}/summary.html?project=scz-demo", 25, False, None, None),
        ]
        for name, url, wait, full, before, sel in jobs:
            if only and name not in only:
                continue
            try:
                shot(page, name, url, wait, full, before, sel)
            except Exception as e:
                print("  FAILED", name, e)
        b.close()


if __name__ == "__main__":
    main(sys.argv[1:] or None)
