"""Interactive screenshots for the intro video (scrolling, clicking): see capture.py."""
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

OUT = Path(__file__).parent / "shots"
DESK = "http://localhost:8765"


def main():
    with sync_playwright() as p:
        b = p.chromium.launch(channel="msedge", headless=True)
        page = b.new_page(viewport={"width": 1920, "height": 1080})

        page.goto(f"{DESK}/viewer.html?project=cad-demo&locus=1", wait_until="domcontentloaded", timeout=180_000)
        time.sleep(14)
        page.get_by_text("LD r²", exact=False).first.scroll_into_view_if_needed()
        page.mouse.wheel(0, 500)
        time.sleep(3)
        page.screenshot(path=str(OUT / "viewer_pcsk9_ld.png"))
        print("  viewer_pcsk9_ld")

        page.goto(f"{DESK}/viewer.html?project=cad-demo&locus=1", wait_until="domcontentloaded", timeout=180_000)
        time.sleep(14)
        page.get_by_text("Show ranking", exact=True).first.click()
        time.sleep(5)
        page.screenshot(path=str(OUT / "viewer_pcsk9_rank.png"))
        print("  viewer_pcsk9_rank")

        page.goto(f"{DESK}/viewer.html?project=t2d-demo&locus=87", wait_until="domcontentloaded", timeout=180_000)
        time.sleep(14)
        page.get_by_text("Show genes", exact=True).first.click()
        time.sleep(4)
        page.screenshot(path=str(OUT / "viewer_tcf7l2_magma.png"))
        print("  viewer_tcf7l2_magma")
        b.close()


if __name__ == "__main__":
    main()
