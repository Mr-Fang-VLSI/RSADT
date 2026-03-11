#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced visualization for multi-column (K=2/4) results.

Now computes total HPWL as:
    Total_HPWL = K * HPWL_eq + (K - 1) * pitch * M

and plots both K=2 and K=4 in one figure:
    - HPWL_eq (solid lines, different colors)
    - Runtime (dashed lines, same color)
"""
import argparse, csv, math, re
from pathlib import Path
from typing import Dict, List
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_csv_rows(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({(k or "").strip(): (v or "").strip() for k,v in r.items()})
    return rows

def to_num(s: str):
    if s is None: return math.nan
    t = s.strip()
    if t=="" or t.lower()=="na": return math.nan
    try:
        if any(c in t for c in (".","e","E")): return float(t)
        return int(t)
    except Exception:
        try: return float(t)
        except Exception: return math.nan

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def plot_overlay_byK(byK, xkey, ykey, title, out_png):
    plt.figure()
    for K in sorted(byK.keys()):
        x = byK[K].get(xkey, []); y = byK[K].get(ykey, [])
        x2,y2=[],[]
        for i,v in enumerate(x):
            if math.isnan(v): continue
            yv = y[i] if i<len(y) else math.nan
            if math.isnan(yv): continue
            x2.append(v); y2.append(yv)
        if x2 and y2:
            plt.plot(x2, y2, label=f"K={K}")
    plt.xlabel(xkey); plt.ylabel(ykey)
    plt.title(title)
    plt.legend(); plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout(); plt.savefig(out_png); plt.close()

def plot_combined_hpwl_runtime(byK, title, out_png):
    """HPWL_eq (solid) and runtime (dashed) for K=2,4 on same plot"""
    plt.figure()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for idx,K in enumerate(sorted(byK.keys())):
        d = byK[K]
        x = d["T"]
        hpwl = d["total_HPWL"]
        runtime = d["total_runtime_sec"]
        xv,yh,yr=[],[],[]
        for i,v in enumerate(x):
            if math.isnan(v): continue
            hv = hpwl[i] if i<len(hpwl) else math.nan
            rv = runtime[i] if i<len(runtime) else math.nan
            if math.isnan(hv) or math.isnan(rv): continue
            xv.append(v); yh.append(hv); yr.append(rv)
        if not xv: continue
        c = colors[idx % len(colors)]
        plt.plot(xv, yh, '-', color=c, label=f"Total_HPWL K={K}")
        plt.plot(xv, yr, '--', color=c, label=f"Runtime K={K}")
    plt.xlabel("T")
    plt.ylabel("Total_HPWL / Runtime(s)")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def process_root(root: Path, pitch: float, out_prefix=""):
    glob = root / "summary_multicol_all.txt"
    if not glob.exists():
        raise SystemExit(f"[ERR] not found: {glob}")
    rows = read_csv_rows(glob)

    byK: Dict[int, Dict[str, List[float]]] = {}
    M = None
    for r in rows:
        try: K = int(r.get("K","0") or "0")
        except: continue
        if K not in (2,4): continue
        if M is None: M = to_num(r.get("M",""))
        d = byK.setdefault(K, {
            "T": [], "total_runtime_sec": [], "last_HPWL_eq": [],
            "tighten_steps": [], "first_meet_T_round": [], "avg_solve_ms": [],
            "total_HPWL": []
        })
        T = to_num(r.get("T",""))
        he = to_num(r.get("last_HPWL_eq",""))
        runtime = to_num(r.get("total_runtime_sec",""))
        first = to_num(r.get("first_meet_T_round",""))
        tight = to_num(r.get("tighten_steps",""))
        avgms = to_num(r.get("avg_solve_ms",""))
        total_hpwl = K * he + (K - 1) * pitch * M
        d["T"].append(T)
        d["last_HPWL_eq"].append(he)
        d["total_runtime_sec"].append(runtime)
        d["first_meet_T_round"].append(first)
        d["tighten_steps"].append(tight)
        d["avg_solve_ms"].append(avgms)
        d["total_HPWL"].append(total_hpwl)

    plots_dir = root / "plots_total"
    ensure_dir(plots_dir)

    # Combined Total HPWL + Runtime overlay
    plot_combined_hpwl_runtime(byK, "Total HPWL (incl cross-column nets) & Runtime vs T (K=2,4)",
                               plots_dir / f"{out_prefix}combined_totalHPWL_runtime_vs_T.png")

    # Other standard overlays
    plot_overlay_byK(byK, "T", "total_HPWL", f"Total HPWL vs T (pitch={pitch})",
                     plots_dir / f"{out_prefix}total_HPWL_vs_T_byK.png")
    plot_overlay_byK(byK, "T", "tighten_steps", "Tighten steps vs T (K=2,4)",
                     plots_dir / f"{out_prefix}global_tighten_steps_vs_T_byK.png")
    plot_overlay_byK(byK, "T", "first_meet_T_round", "First meet T vs T (K=2,4)",
                     plots_dir / f"{out_prefix}global_first_meet_round_vs_T_byK.png")

    print(f"[OK] Plots saved under: {plots_dir}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result_root", type=Path, help="Path to results_multicol_<MxH>_<timestamp>")
    ap.add_argument("--pitch", type=float, default=15.0, help="Column pitch (default=15)")
    ap.add_argument("--out-prefix", default="", help="Prefix for filenames")
    args = ap.parse_args()
    if not args.result_root.exists():
        raise SystemExit(f"[ERR] not found: {args.result_root}")
    process_root(args.result_root, args.pitch, args.out_prefix)

if __name__ == "__main__":
    main()
