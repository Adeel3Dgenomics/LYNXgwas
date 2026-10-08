"""Store the GRCh38 coordinates of compact scE2G files as offsets from GRCh37 (23% smaller, nothing lost).

    d38    = start38 - start
    dtss38 = (tss38 - tss) - d38
    dlen38 = (end38 - start38) - (end - start)
so start38 = start + d38, end38 = start38 + (end - start) + dlen38, tss38 = tss + d38 + dtss38.
Mostly zeros and one repeated offset per region, which gzip compresses far better than absolute positions.

    python scripts/regulatory/delta38.py regulatory_links/sce2g        (converts in place; updates manifest.json)
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

DIR = Path(sys.argv[1])
COLS = ["chr", "start", "end", "class", "gene", "tss", "self_promoter", "e2g", "e2g_raw", "abc", "kendall", "arc",
        "d38", "dtss38", "dlen38"]


def convert(path: Path) -> None:
    text = gzip.open(path, "rt").read().split("\n")
    hdr = [ln for ln in text if ln.startswith("#")]
    rows = [ln.split("\t") for ln in text if ln and not ln.startswith("#")]
    h = rows[0]
    if "d38" in h:
        return
    i = {c: h.index(c) for c in h}
    out = [ln.replace("GRCh37 (lifted) + GRCh38", "GRCh37 (lifted) + GRCh38 as offsets: start38 = start + d38, "
                      "end38 = start38 + (end - start) + dlen38, tss38 = tss + d38 + dtss38") for ln in hdr]
    out.append("\t".join(COLS))
    for r in rows[1:]:
        s, e, t = int(r[i["start"]]), int(r[i["end"]]), int(r[i["tss"]])
        s8, e8, t8 = int(r[i["start38"]]), int(r[i["end38"]]), int(r[i["tss38"]])
        d = s8 - s
        out.append("\t".join([r[i[c]] for c in COLS[:12]] + [str(d), str(t8 - t - d), str((e8 - s8) - (e - s))]))
    path.write_bytes(gzip.compress(("\n".join(out) + "\n").encode("utf-8"), 9))


def main():
    mf = DIR / "manifest.json"
    m = json.loads(mf.read_text())
    for ct in m["cell_types"]:
        p = DIR / ct["file"]
        convert(p)
        ct["bytes"] = p.stat().st_size
        ct["md5"] = hashlib.md5(p.read_bytes()).hexdigest()
    m["columns"] = COLS
    m["builds"] = ("GRCh37 (lifted with the UCSC hg38ToHg19 chain) and the original GRCh38, stored as offsets "
                   "(start38 = start + d38, end38 = start38 + (end - start) + dlen38, tss38 = tss + d38 + dtss38)")
    mf.write_text(json.dumps(m, indent=1))
    print(f"{len(m['cell_types'])} files, {sum(c['bytes'] for c in m['cell_types']) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
