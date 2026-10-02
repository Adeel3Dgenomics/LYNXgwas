"""Extract an LD sub-matrix for given SNPs from the UK Biobank LD matrices (Weissbrod et al. 2020,
https://registry.opendata.aws/ukbb-ld/; 337k British-ancestry individuals, GRCh37, 3 Mb windows).

A 1000 Genomes EUR panel has 503 people; for GWAS with hundreds of thousands of samples its LD
errors look like extra signals to SuSiE (all L effects used, COJO runaway). The UKB matrices are a
much closer stand-in for large European GWAS.

Usage:
  python ukbb_ld.py <ukbb_dir> <snps.tsv> <out_prefix>
    snps.tsv   : chr<TAB>pos<TAB>a1<TAB>a2 per SNP (a1 = allele the z-scores refer to), with header
    out_prefix : writes <out_prefix>.ld (space-separated matrix of the matched SNPs, in input order)
                 and <out_prefix>.idx (1-based input row of each matched SNP)
Exit code 2 when no window covers the SNPs (caller falls back to the reference panel).

Windows are <ukbb_dir>/chr{C}_{S}_{S+3000000}.npz (+ .gz SNP list), S = 1, 1000001, ...
Sign convention: correlations in a window share one allele coding, so flipping both SNPs of a
pair leaves r unchanged; a SNP whose (a1, a2) is swapped relative to the window flips sign.
"""
import gzip, os, sys

import numpy as np
import scipy.sparse as sp

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def pick_window(ukbb_dir, chrom, lo, hi):
    """The available window containing [lo, hi] whose centre is closest to the SNPs' centre."""
    best = None
    for f in os.listdir(ukbb_dir):
        if not f.endswith(".npz") or not f.startswith(f"chr{chrom}_"):
            continue
        _, s, e = f[:-4].split("_")
        s, e = int(s), int(e)
        if s <= lo and hi < e and os.path.isfile(os.path.join(ukbb_dir, f[:-4] + ".gz")):
            d = abs((s + e) / 2 - (lo + hi) / 2)
            if best is None or d < best[0]:
                best = (d, f[:-4])
    return None if best is None else os.path.join(ukbb_dir, best[1])


def main(ukbb_dir, snp_file, out_prefix):
    rows = []
    with open(snp_file) as f:
        next(f)
        for line in f:
            c, p, a1, a2 = line.split()[:4]
            rows.append((c.replace("chr", ""), int(p), a1.upper(), a2.upper()))
    if not rows:
        sys.exit(2)
    chrom = rows[0][0]
    win = pick_window(ukbb_dir, chrom, min(r[1] for r in rows), max(r[1] for r in rows))
    if win is None:
        print("no UKB LD window covers these SNPs", file=sys.stderr)
        sys.exit(2)

    # Window SNP list: rsid chromosome position allele1 allele2
    pos_index = {}
    with gzip.open(win + ".gz", "rt") as f:
        header = f.readline().split()
        ip, i1, i2 = header.index("position"), header.index("allele1"), header.index("allele2")
        for k, line in enumerate(f):
            x = line.split()
            pos_index.setdefault(int(x[ip]), []).append((k, x[i1].upper(), x[i2].upper()))

    matched, sign, idx = [], [], []
    for i, (c, p, a1, a2) in enumerate(rows, start=1):
        for k, u1, u2 in pos_index.get(p, []):
            if (a1, a2) == (u1, u2) or (COMP.get(a1), COMP.get(a2)) == (u1, u2):
                matched.append(k); sign.append(1.0); idx.append(i); break
            if (a1, a2) == (u2, u1) or (COMP.get(a1), COMP.get(a2)) == (u2, u1):
                matched.append(k); sign.append(-1.0); idx.append(i); break

    M = sp.load_npz(win + ".npz").tocsr()
    sub = M[matched][:, matched].toarray()
    R = sub + sub.T                       # stored as one triangle
    np.fill_diagonal(R, 1.0)
    s = np.array(sign)
    R = R * np.outer(s, s)
    np.savetxt(out_prefix + ".ld", R, fmt="%.6f")
    with open(out_prefix + ".idx", "w") as f:
        f.write("\n".join(map(str, idx)) + "\n")
    print(f"UKB LD: {os.path.basename(win)}, matched {len(idx)}/{len(rows)} SNPs")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
