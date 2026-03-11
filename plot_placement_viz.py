#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_placement_viz.py
=====================

Visualize placement solutions (e.g. placement16_Col_4_T_12.txt).

For each file, draw two figures:
1) Physical layout: DSP site grid, label each cell with logical (i,j).
2) Logical layout: logical array grid, label each cell with its physical y.

Usage:
  python3 plot_placement_viz.py <placement1.txt> [placement2.txt ...]
  python3 plot_placement_viz.py <directory>

Output:
  viz_placements/<file_stem>/
    ├─ *_phys.png   (physical layout)
    └─ *_logic.png  (logical layout)

Dependencies: matplotlib
"""
import sys, re, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Regex definitions ---
HEADER_RE = re.compile(r"(\d+)\s*x\s*(\d+).*?on\s+(\d+)\s+DSP", re.IGNORECASE)
LINE_RE = re.compile(r"([XY])_(\d+)_(\d+)\s+(-?\d+)")

# --- Helpers ---
def parse_file(fp: Path):
    """Return m,h,K,map[(i,j)]=(x,y),max_y"""
    text = fp.read_text(encoding="utf-8", errors="ignore").splitlines()
    m = h = K = None
    data = {}
    for line in text:
        m0 = HEADER_RE.search(line)
        if m0 and (m is None or h is None or K is None):
            m, h, K = map(int, m0.groups())
        m1 = LINE_RE.match(line.strip())
        if m1:
            kind, i, j, v = m1.groups()
            i, j, v = int(i), int(j), int(v)
            data.setdefault((i, j), {})[kind] = v
    coords = {}
    max_y = 0
    for (i, j), d in data.items():
        if "X" in d and "Y" in d:
            coords[(i, j)] = (d["X"], d["Y"])
            max_y = max(max_y, d["Y"])
    if m is None: m = max(i for i, _ in coords) + 1
    if h is None: h = max(j for _, j in coords) + 1
    if K is None: K = max(x for x, _ in coords.values()) + 1
    return m, h, K, coords, max_y

def draw_grid(ax, w, h):
    for x in range(w+1):
        ax.plot([x, x], [0, h], lw=0.8)
    for y in range(h+1):
        ax.plot([0, w], [y, y], lw=0.8)
    ax.set_xlim(0, w); ax.set_ylim(0, h)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.axis("off")

def annotate(ax, x, y, text, size=6):
    ax.text(x, y, text, ha="center", va="center", fontsize=size)

def fig_size(w, h, base=0.4):
    return max(3, w*base), min(30, max(3, h*base))

# --- Plotters ---
def plot_physical(fp, m, h, K, coords, max_y, outdir):
    """Physical view: grid=(K, max_y+1), label=(i,j)."""
    H = max_y + 1
    fig, ax = plt.subplots(figsize=fig_size(K, H))
    draw_grid(ax, K, H)
    for (i,j),(x,y) in coords.items():
        if 0<=x<K and 0<=y<H:
            annotate(ax, x+0.5, y+0.5, f"({i},{j})")
    ax.set_title(f"Physical layout\n{fp.name} (K={K}, height={H})")
    out = outdir / f"{fp.stem}_phys.png"
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)
    print(f"[OK] Physical -> {out}")

def plot_logical(fp, m, h, K, coords, outdir):
    """Logical view: grid=(h,m), label=y."""
    fig, ax = plt.subplots(figsize=fig_size(h, m))
    draw_grid(ax, h, m)
    for i in range(m):
        for j in range(h):
            y = coords.get((i,j), (None,None))[1]
            if y is not None:
                annotate(ax, j+0.5, i+0.5, str(y))
    ax.set_title(f"Logical layout (labels=y)\n{fp.name}")
    out = outdir / f"{fp.stem}_logic.png"
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)
    print(f"[OK] Logical -> {out}")

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(2)
    inputs = []
    for a in sys.argv[1:]:
        p = Path(a)
        if p.is_dir():
            inputs += sorted(p.glob("placement*.txt"))
        elif p.is_file():
            inputs.append(p)
    if not inputs:
        print("[ERR] No placement files found."); sys.exit(1)
    outroot = Path("viz_placements"); outroot.mkdir(exist_ok=True)
    for fp in inputs:
        try:
            m,h,K,coords,max_y = parse_file(fp)
            outdir = outroot / fp.stem; outdir.mkdir(exist_ok=True)
            plot_physical(fp,m,h,K,coords,max_y,outdir)
            plot_logical(fp,m,h,K,coords,outdir)
        except Exception as e:
            print(f"[ERR] {fp}: {e}")

if __name__ == "__main__":
    main()
