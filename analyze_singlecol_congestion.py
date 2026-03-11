#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_singlecol_congestion.py
--------------------------------
For each placement*.txt file:
  1) Parse mapping (i,j) -> (x,y)
  2) Select one physical column (default x=0)
  3) Build intra-strip neighbor nets: vertical (i+1,j) and horizontal (i,j+1) if both endpoints in the same column
  4) Plot:
     - Histogram of net lengths: |Δy|
     - Congestion profile along y: for each segment [y, y+1), how many nets cover it (interval overlap count)

Usage:
  python3 analyze_singlecol_congestion.py <file_or_dir> [more ...] --col 0 --bins 50 --dpi 180

Notes:
  * No seaborn. Each chart has its own figure. No custom colors (per your plotting rules).
  * Cross-seam nets (between different physical columns) are ignored here; we analyze a single strip.
"""

import argparse
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HEADER_RE = re.compile(r"(\d+)\s*x\s*(\d+).*?on\s+(\d+)\s+DSP", re.IGNORECASE)
LINE_RE   = re.compile(r"([XY])_(\d+)_(\d+)\s+(-?\d+)")

def parse_placement(fp: Path):
    """Return dict: coords[(i,j)] = (x,y); and inferred (m,h,K,max_y) if available."""
    text = fp.read_text(encoding="utf-8", errors="ignore").splitlines()
    coords = {}
    m = h = K = None
    max_y = 0
    for line in text:
        m0 = HEADER_RE.search(line)
        if m0 and (m is None or h is None or K is None):
            m, h, K = map(int, m0.groups())
        m1 = LINE_RE.match(line.strip())
        if m1:
            kind, i, j, v = m1.groups()
            i, j, v = int(i), int(j), int(v)
            if (i, j) not in coords:
                coords[(i, j)] = [None, None]  # x, y
            if kind == "X":
                coords[(i, j)][0] = v
            else:
                coords[(i, j)][1] = v
            if kind == "Y" and v > max_y:
                max_y = v
    # filter complete entries
    coords = {k: (xy[0], xy[1]) for k, xy in coords.items() if xy[0] is not None and xy[1] is not None}
    if m is None: m = max(i for i, _ in coords.keys()) + 1
    if h is None: h = max(j for _, j in coords.keys()) + 1
    if K is None: K = max(x for x, _ in coords.values()) + 1
    return coords, m, h, K, max_y

def collect_strip(coords, x_col):
    """Return nodes and neighbor nets (within the same physical column x_col)."""
    # nodes in this column
    nodes = {ij for ij, (x, y) in coords.items() if x == x_col}
    # build adjacency nets (i+1,j) and (i,j+1) if both in this column
    nets = []
    for (i, j) in nodes:
        if (i+1, j) in nodes:
            nets.append(((i, j), (i+1, j)))
        if (i, j+1) in nodes:
            nets.append(((i, j), (i, j+1)))
    return nodes, nets

def net_lengths(coords, nets):
    """Return list of |Δy| for given nets."""
    lens = []
    for (u, v) in nets:
        yu = coords[u][1]; yv = coords[v][1]
        lens.append(abs(yu - yv))
    return lens

def coverage_profile(coords, nets):
    """Return coverage per unit segment [y, y+1)."""
    # Determine height from max y among endpoints
    if not coords:
        return []
    H = max(y for (_, y) in coords.values()) + 1
    if H <= 1:
        return [0]
    diff = [0]*(H)  # difference array; coverage on segments uses 0..H-2
    for (u, v) in nets:
        yu = coords[u][1]; yv = coords[v][1]
        a, b = (yu, yv) if yu <= yv else (yv, yu)
        if a == b:
            continue
        # cover segments [a, b): add 1
        diff[a] += 1
        diff[b] -= 1
    cov = [0]*(H-1)
    run = 0
    for y in range(H-1):
        run += diff[y]
        cov[y] = run
    return cov

def plot_histogram(lens, out_png, title):
    if len(lens) == 0:
        # still create an empty plot for consistency
        fig, ax = plt.subplots(figsize=(6,4))
        ax.set_title(title + " (no nets)")
        ax.set_xlabel("|Δy| (net length)")
        ax.set_ylabel("#nets")
        fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)
        return
    max_len = max(lens)
    bins = min(80, max(10, max_len))  # auto bins
    fig, ax = plt.subplots(figsize=(7,4.5))
    ax.hist(lens, bins=bins)
    ax.set_title(title)
    ax.set_xlabel("|Δy| (net length along the column)")
    ax.set_ylabel("# nets")
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)

def plot_coverage(cov, out_png, title):
    fig, ax = plt.subplots(figsize=(8,3.5))
    ax.plot(range(len(cov)), cov)
    ax.set_title(title)
    ax.set_xlabel("segment index (between y and y+1)")
    ax.set_ylabel("# nets covering segment")
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)

def process_file(fp: Path, x_col: int, out_root: Path):
    coords, m, h, K, max_y = parse_placement(fp)
    if not coords:
        print(f"[WARN] {fp}: no (X,Y) pairs parsed.")
        return
    nodes, nets = collect_strip(coords, x_col)
    # filter to nets whose both endpoints exist in coords (already ensured)
    lens = net_lengths(coords, nets)
    cov  = coverage_profile(coords, nets)

    outdir = out_root / f"{fp.stem}_col{x_col}"
    outdir.mkdir(parents=True, exist_ok=True)

    # Plots
    plot_histogram(lens, outdir / f"{fp.stem}_col{x_col}_netlen_hist.png",
                   f"Net length histogram (x={x_col})\n{fp.name}")
    plot_coverage(cov, outdir / f"{fp.stem}_col{x_col}_coverage.png",
                  f"Segment coverage along column (x={x_col})\n{fp.name}")

    # Console summary
    if lens:
        lens_sorted = sorted(lens)
        n = len(lens_sorted)
        p50 = lens_sorted[n//2]
        p95 = lens_sorted[int(0.95*(n-1))]
        print(f"[OK] {fp.name}  x={x_col}  nets={n}  len[max]={max(lens_sorted)}  p50={p50}  p95={p95}  "
              f"coverage[max]={max(cov) if cov else 0}")
    else:
        print(f"[OK] {fp.name}  x={x_col}  nets=0")

def collect_inputs(paths):
    files = []
    for p in paths:
        P = Path(p)
        if P.is_dir():
            files.extend(sorted(P.glob("placement*.txt")))
        elif P.is_file():
            files.append(P)
    return files

def main():
    ap = argparse.ArgumentParser(description="Analyze single-column net length & congestion from placement*.txt")
    ap.add_argument("inputs", nargs="+", help="placement file(s) or directory(ies)")
    ap.add_argument("--col", type=int, default=0, help="physical column x to analyze (default 0)")
    ap.add_argument("--out", type=str, default="viz_congestion_singlecol", help="output root dir")
    args = ap.parse_args()

    files = collect_inputs(args.inputs)
    if not files:
        print("[ERR] No placement files found.")
        return
    out_root = Path(args.out); out_root.mkdir(exist_ok=True)

    for fp in files:
        process_file(fp, args.col, out_root)

if __name__ == "__main__":
    main()
