"""Compact one scE2G prediction file into the LYNXgwas regulatory-links format.

scE2G (Sheth, Qiu et al., Engreitz lab; https://data.broadinstitute.org/scE2Gpredictions/) predicts, per cell
type, every candidate element-gene pair genome-wide (GRCh38, 27 columns, ~700 MB per cell type). LYNXgwas needs
only the pairs scE2G itself calls (E2G.Score.qnorm >= 0.177, the scE2G_multiome threshold) and a few columns:

    chr  start  end  class  gene  tss  self_promoter  e2g  e2g_raw  abc  kendall  arc  start38  end38  tss38

`start/end/tss` are GRCh37 (lifted with the UCSC hg38ToHg19 chain; pairs whose element or TSS does not lift,
or lifts to another chromosome or to a span more than 2x its length, are dropped and counted), `*38` the original
GRCh38 positions. Scores are rounded to 4 significant digits. Output: one gzipped TSV per cell type, sorted by
chr and start, with a header line.

    python scripts/regulatory/compact_sce2g.py <raw.tsv.gz> <cell type label> <hg38ToHg19.over.chain.gz> <out.tsv.gz>
"""
import csv
import gzip
import sys
import time

from pyliftover import LiftOver

THR = 0.177
KEEP_OUT = ["chr", "start", "end", "class", "gene", "tss", "self_promoter", "e2g", "e2g_raw", "abc", "kendall", "arc",
            "start38", "end38", "tss38"]


def sig(x: str, n: int = 4) -> str:
    try:
        v = float(x)
    except ValueError:
        return ""
    return f"{v:.{n}g}"


def lift(lo, chrom, pos):
    r = lo.convert_coordinate(chrom, pos)
    if not r or r[0][0] != chrom:
        return None
    return r[0][1]


def main(raw, label, chain, out):
    t0 = time.time()
    lo = LiftOver(chain)
    rows, n_in, n_kept, n_unlifted = [], 0, 0, 0
    cache = {}
    with gzip.open(raw, "rt", newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        h = next(rd)
        i = {c: h.index(c) for c in ["chr", "start", "end", "class", "TargetGene", "TargetGeneTSS", "isSelfPromoter",
                                     "E2G.Score.qnorm", "E2G.Score", "ABC.Score", "Kendall", "ARC.E2G.Score"]}
        for r in rd:
            n_in += 1
            try:
                q = float(r[i["E2G.Score.qnorm"]])
            except ValueError:
                continue
            if q < THR:
                continue
            n_kept += 1
            ch, s38, e38, tss38 = r[i["chr"]], int(r[i["start"]]), int(r[i["end"]]), int(float(r[i["TargetGeneTSS"]]))
            key = (ch, s38, e38)
            if key not in cache:
                s37, e37 = lift(lo, ch, s38), lift(lo, ch, e38 - 1)
                ok = s37 is not None and e37 is not None and 0 < (e37 + 1 - s37) <= 2 * (e38 - s38)
                cache[key] = (s37, e37 + 1) if ok else None
            el = cache[key]
            t37 = lift(lo, ch, tss38)
            if el is None or t37 is None:
                n_unlifted += 1
                continue
            rows.append([ch, el[0], el[1], r[i["class"]], r[i["TargetGene"]], t37,
                         "1" if r[i["isSelfPromoter"]] == "True" else "0", sig(r[i["E2G.Score.qnorm"]]),
                         sig(r[i["E2G.Score"]]), sig(r[i["ABC.Score"]]), sig(r[i["Kendall"]]), sig(r[i["ARC.E2G.Score"]]),
                         s38, e38, tss38])
    chrom_order = {f"chr{c}": k for k, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
    rows.sort(key=lambda x: (chrom_order.get(x[0], 99), x[1], x[4]))
    with gzip.open(out, "wt", newline="", compresslevel=9) as g:
        g.write(f"# scE2G_multiome {label}; Sheth, Qiu et al. (Engreitz lab), data.broadinstitute.org/scE2Gpredictions; "
                f"E2G.Score.qnorm >= {THR}; GRCh37 (lifted) + GRCh38\n")
        w = csv.writer(g, delimiter="\t", lineterminator="\n")
        w.writerow(KEEP_OUT)
        w.writerows(rows)
    print(f"{label}: {n_in:,} pairs read, {n_kept:,} >= {THR}, {n_unlifted:,} not lifted, {len(rows):,} written "
          f"({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:5])
