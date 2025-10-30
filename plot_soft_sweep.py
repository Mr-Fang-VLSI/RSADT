#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process results produced by run_soft_sweep.sh (soft policy only).

Given a RESULT_ROOT directory (e.g., results_soft_2025xxxx_xxxxxx), this script:
  1) For each case folder (soft_<M>x<H>_T<Tmin>-<Tmax>): 
     - Loads summary_soft_<M>x<H>.txt
     - Generates overview CSV (overview_<M>x<H>.csv) for T vs key metrics
     - Plots:
       * T_vs_runtime.png  (T on X, total_runtime_sec on Y)
       * T_vs_HPWL_eq.png  (T on X, last_HPWL_eq on Y)
       * Tprime_vs_round_allT.png  (overlay T_prime vs round for all T runs)
       * maxL_vs_round_allT.png    (overlay maxL vs round for all T runs)
       * first_meet_round_vs_T.png (first_meet_T_round vs T)
       * tighten_steps_vs_T.png    (tighten_steps vs T)
       * avg_solve_ms_vs_T.png     (average per-round solve_ms vs T)
       * HPWL_ratio_vs_T.png       (HPWL_w / HPWL_eq vs T)
  2) If pandas is available, also writes an Excel workbook collecting per-case overviews.

Usage:
    python3 plot_soft_sweep.py <RESULT_ROOT> [--out-prefix PFX]

Notes:
- This script uses only matplotlib (no seaborn) and default styles/colors.
- Each chart is created as a single figure (no subplots), per requirement.
"""
import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, List, Any, Tuple

# Optional pandas for Excel export
try:
    import pandas as pd
    HAVE_PANDAS = True
except Exception:
    HAVE_PANDAS = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({k.strip(): v.strip() for k, v in r.items()})
    return rows


def to_num(s: str):
    if s is None:
        return math.nan
    t = s.strip()
    if t == "" or t.lower() == "na":
        return math.nan
    try:
        if any(c in t for c in (".", "e", "E")):
            return float(t)
        return int(t)
    except Exception:
        try:
            return float(t)
        except Exception:
            return math.nan


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def find_case_dirs(root: Path) -> List[Path]:
    # case folders created by run_soft_sweep.sh: soft_<M>x<H>_T<Tmin>-<Tmax>
    return sorted([p for p in root.iterdir() if p.is_dir() and p.name.startswith("soft_")])


def parse_MH_from_case(case_dir: Path) -> Tuple[int, int]:
    m = h = None
    m = re.search(r"soft_(\d+)x(\d+)", case_dir.name)
    if m:
        return int(m.group(1)), int(m.group(2))
    return -1, -1


def plot_xy(x, y, xlabel: str, ylabel: str, title: str, out_png: Path):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def parse_tprime_file(fp: Path):
    """
    Parse tprime_evolution_soft_<MxH>_T_<T>.txt (round and tighten entries).
    Returns:
      T (int), rounds dict with arrays for keys: round, T_prime, maxL, WNS_T, HPWL_eq, HPWL_w, solve_ms
    """
    T_match = re.search(r"_T_(\d+)", fp.name)
    T_val = int(T_match.group(1)) if T_match else None

    rounds = {"round": [], "T_prime": [], "maxL": [], "WNS_T": [], "HPWL_eq": [], "HPWL_w": [], "solve_ms": []}
    with fp.open("r", encoding="utf-8") as f:
        # skip comment lines
        reader = csv.DictReader((line for line in f if not line.startswith("#")))
        for row in reader:
            if row.get("type", "") != "round":
                continue
            rounds["round"].append(to_num(row.get("round", "")))
            rounds["T_prime"].append(to_num(row.get("T_prime", "")))
            rounds["maxL"].append(to_num(row.get("maxL", "")))
            rounds["WNS_T"].append(to_num(row.get("WNS_T", "")))
            rounds["HPWL_eq"].append(to_num(row.get("HPWL_eq", "")))
            rounds["HPWL_w"].append(to_num(row.get("HPWL_w", "")))
            rounds["solve_ms"].append(to_num(row.get("solve_ms", "")))
    return T_val, rounds


def plot_overlay_by_T(series_by_T: Dict[int, Dict[str, List[float]]], ykey: str, title: str, out_png: Path):
    """
    Overlay multiple T curves (one per requested T) for ykey vs round.
    """
    plt.figure()
    for T, d in sorted(series_by_T.items()):
        rr = [v for v in d["round"] if not math.isnan(v)]
        yy_raw = d.get(ykey, [])
        yy = [yy_raw[i] for i, v in enumerate(d["round"]) if not math.isnan(v)]
        if rr and yy:
            plt.plot(rr, yy, label=f"T={T}")
    plt.xlabel("round")
    plt.ylabel(ykey)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def process_case(case_dir: Path, out_prefix: str = ""):
    m, h = parse_MH_from_case(case_dir)
    if m <= 0:
        print(f"[WARN] Skip unknown case dir: {case_dir}")
        return

    plots_dir = case_dir / "plots"
    ensure_dir(plots_dir)

    # 1) Load per-case summary
    summ_file = next(case_dir.glob("summary_soft_*x*.txt"), None)
    if not summ_file or not summ_file.exists():
        print(f"[WARN] summary not found in {case_dir}")
        return
    rows = read_csv_rows(summ_file)

    # Build overview table (ordered by T)
    overview_cols = [
        "M","H","T","rounds_total","total_runtime_sec",
        "last_round","last_maxL","last_HPWL_eq","last_HPWL_w","last_WNS_T","meet_T","first_meet_T_round",
        "init_Tprime","tighten_steps","final_Tprime","min_maxL","round_of_min","max_maxL","round_of_max",
        "avg_solve_ms","last_solve_ms","last_picked",
        "alpha","eta","gamma","beta","sigma","shrink","cap","scale","dV","pow2","freeze",
        "log_path","placement_path","tprime_path"
    ]
    # Some rows might have "NA" – handle conversion later in plotting
    def keyT(r): 
        try: return int(r.get("T","0"))
        except: return 0
    rows_sorted = sorted(rows, key=keyT)

    # write overview CSV
    overview_path = case_dir / f"{out_prefix}overview_{m}x{h}.csv"
    with overview_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=overview_cols)
        w.writeheader()
        for r in rows_sorted:
            w.writerow({k: r.get(k, "") for k in overview_cols})

    # 2) Plots from overview
    T_vals = [to_num(r.get("T","")) for r in rows_sorted]
    runtime = [to_num(r.get("total_runtime_sec","")) for r in rows_sorted]
    hpwl_eq = [to_num(r.get("last_HPWL_eq","")) for r in rows_sorted]
    first_meet = [to_num(r.get("first_meet_T_round","")) for r in rows_sorted]
    tighten_steps = [to_num(r.get("tighten_steps","")) for r in rows_sorted]
    avg_ms = [to_num(r.get("avg_solve_ms","")) for r in rows_sorted]
    hpwl_ratio = []
    for r in rows_sorted:
        he = to_num(r.get("last_HPWL_eq",""))
        hw = to_num(r.get("last_HPWL_w",""))
        hpwl_ratio.append(hw / he if (not math.isnan(he) and he != 0 and not math.isnan(hw)) else math.nan)

    plot_xy(T_vals, runtime, "T", "total_runtime_sec", f"Runtime vs T ({m}x{h})", plots_dir / f"{out_prefix}T_vs_runtime.png")
    plot_xy(T_vals, hpwl_eq, "T", "last_HPWL_eq", f"HPWL_eq vs T ({m}x{h})", plots_dir / f"{out_prefix}T_vs_HPWL_eq.png")
    plot_xy(T_vals, first_meet, "T", "first_meet_T_round", f"First meet T vs T ({m}x{h})", plots_dir / f"{out_prefix}first_meet_round_vs_T.png")
    plot_xy(T_vals, tighten_steps, "T", "tighten_steps", f"Tighten steps vs T ({m}x{h})", plots_dir / f"{out_prefix}tighten_steps_vs_T.png")
    plot_xy(T_vals, avg_ms, "T", "avg_solve_ms", f"Avg DP solve_ms vs T ({m}x{h})", plots_dir / f"{out_prefix}avg_solve_ms_vs_T.png")
    plot_xy(T_vals, hpwl_ratio, "T", "HPWL_w/HPWL_eq", f"HPWL ratio vs T ({m}x{h})", plots_dir / f"{out_prefix}HPWL_ratio_vs_T.png")

    # 3) Overlay T_prime vs round and maxL vs round across all T runs
    t_files = sorted(case_dir.glob("tprime_evolution_soft_*_T_*.txt"))
    series_by_T = {}
    for fp in t_files:
        T_val, series = parse_tprime_file(fp)
        if T_val is not None:
            series_by_T[T_val] = series

    if series_by_T:
        plot_overlay_by_T(series_by_T, "T_prime", f"T_prime vs round ({m}x{h})", plots_dir / f"{out_prefix}Tprime_vs_round_allT.png")
        plot_overlay_by_T(series_by_T, "maxL", f"maxL vs round ({m}x{h})", plots_dir / f"{out_prefix}maxL_vs_round_allT.png")

    print(f"[OK] Case {m}x{h}: wrote {overview_path}, plots in {plots_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result_root", type=Path, help="Path to results_soft_YYYYmmdd_HHMMSS (output of run_soft_sweep.sh)")
    ap.add_argument("--out-prefix", default="", help="Optional prefix for output filenames")
    args = ap.parse_args()

    root: Path = args.result_root
    if not root.exists():
        raise SystemExit(f"[ERR] result_root not found: {root}")

    # Either a root with case subdirs, or a single case dir
    cases = [p for p in root.iterdir() if p.is_dir() and p.name.startswith("soft_")] if root.is_dir() else []
    if not cases and root.name.startswith("soft_"):
        cases = [root]

    if not cases:
        raise SystemExit(f"[ERR] No case dirs (soft_*) found in {root}")

    for case_dir in cases:
        process_case(case_dir, args.out_prefix)

    # Optional: aggregate Excel workbook
    if HAVE_PANDAS and len(cases) > 0:
        xlsx = root / (args.out_prefix + "soft_overview.xlsx")
        with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
            for case_dir in cases:
                m, h = parse_MH_from_case(case_dir)
                ov = case_dir / f"{args.out_prefix}overview_{m}x{h}.csv"
                if ov.exists():
                    df = pd.read_csv(ov)
                    df.to_excel(writer, sheet_name=f"{m}x{h}", index=False)
        print(f"[OK] Wrote Excel: {xlsx}")


if __name__ == "__main__":
    main()
