#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot T_prime evolution and related metrics for each policy from
`tprime_evolution_<policy>.txt` files produced by run_compare_any.sh.

Usage:
    python3 plot_tprime.py <RESULT_DIR> [--out-prefix PFX]

This will:
  1) Parse each tprime_evolution_<policy>.txt under RESULT_DIR
  2) Save three figures into RESULT_DIR:
     - tprime_vs_round.png
     - maxL_vs_round.png
     - wnsT_vs_round.png
  3) Save a CSV summary into RESULT_DIR:
     - tprime_summary.csv
"""
import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple, Any

# Optional: use pandas if available; otherwise fallback to stdlib csv
try:
    import pandas as pd
    HAVE_PANDAS = True
except Exception:
    HAVE_PANDAS = False

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt


def _to_num(x: str):
    """Convert string to float (or int if integer-looking). Return NaN on empty."""
    if x is None:
        return math.nan
    s = str(x).strip()
    if s == "" or s.lower() == "na":
        return math.nan
    try:
        if "." in s or "e" in s.lower():
            return float(s)
        return int(s)
    except Exception:
        try:
            return float(s)
        except Exception:
            return math.nan


def parse_tprime_file(path: Path) -> Tuple[dict, dict, dict]:
    """Parse a tprime_evolution_<policy>.txt file.
    Returns:
        meta: dict (unused placeholder)
        rounds: dict of lists for type=='round'
        tighten: dict of lists for type=='tighten'
    """
    meta = {}
    rounds = {
        "round": [], "T_prime": [], "T_target": [], "maxL": [], "WNS_Tprime": [], "WNS_T": [],
        "picked": [], "HPWL_eq": [], "HPWL_w": [], "solve_ms": []
    }
    tighten = {"round": [], "old_T_prime": [], "new_T_prime": [], "sigma_maxL": [], "shrink_val": []}

    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader((line for line in f if not line.startswith("#")))
        header = None
        for row in reader:
            if not row:
                continue
            if header is None:
                header = row
                continue
            rowd = {header[i]: row[i] if i < len(row) else "" for i in range(len(header))}
            rtype = rowd.get("type", "").strip()

            if rtype == "round":
                rounds["round"].append(_to_num(rowd.get("round", "")))
                rounds["T_prime"].append(_to_num(rowd.get("T_prime", "")))
                rounds["T_target"].append(_to_num(rowd.get("T_target", "")))
                rounds["maxL"].append(_to_num(rowd.get("maxL", "")))
                rounds["WNS_Tprime"].append(_to_num(rowd.get("WNS_Tprime", "")))
                rounds["WNS_T"].append(_to_num(rowd.get("WNS_T", "")))
                rounds["picked"].append(_to_num(rowd.get("picked", "")))
                rounds["HPWL_eq"].append(_to_num(rowd.get("HPWL_eq", "")))
                rounds["HPWL_w"].append(_to_num(rowd.get("HPWL_w", "")))
                rounds["solve_ms"].append(_to_num(rowd.get("solve_ms", "")))
            elif rtype == "tighten":
                tighten["round"].append(_to_num(rowd.get("round", "")))
                tighten["old_T_prime"].append(_to_num(rowd.get("old_T_prime", "")))
                tighten["new_T_prime"].append(_to_num(rowd.get("new_T_prime", "")))
                tighten["sigma_maxL"].append(_to_num(rowd.get("sigma_maxL", "")))
                tighten["shrink_val"].append(_to_num(rowd.get("shrink_val", "")))

    return meta, rounds, tighten


def plot_series_by_policy(data: Dict[str, dict], ykey: str, out_png: Path, title: str):
    """Overlay ykey vs round for each policy in one figure (no explicit colors)."""
    plt.figure()
    for pol, d in sorted(data.items()):
        if not d["round"]:
            continue
        xs = [v for v in d["round"] if not math.isnan(v)]
        ys_raw = d.get(ykey, [])
        ys = [ys_raw[i] for i, v in enumerate(d["round"]) if not math.isnan(v)]
        if xs and ys:
            plt.plot(xs, ys, label=pol)
    plt.xlabel("round")
    plt.ylabel(ykey)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def add_tighten_vlines(ax, tighten_rounds: List[float]):
    for r in tighten_rounds:
        try:
            ax.axvline(r, linestyle=":", linewidth=0.8)
        except Exception:
            pass


def plot_with_tighten_marks(data: Dict[str, dict], tighten: Dict[str, dict], ykey: str, out_png: Path, title: str):
    fig = plt.figure()
    ax = fig.gca()
    for pol, d in sorted(data.items()):
        xs = [v for v in d["round"] if not math.isnan(v)]
        ys_raw = d.get(ykey, [])
        ys = [ys_raw[i] for i, v in enumerate(d["round"]) if not math.isnan(v)]
        if xs and ys:
            ax.plot(xs, ys, label=pol)
    all_tighten = []
    for pol, td in tighten.items():
        all_tighten.extend([r for r in td.get("round", []) if not math.isnan(r)])
    add_tighten_vlines(ax, all_tighten)
    ax.set_xlabel("round")
    ax.set_ylabel(ykey)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result_dir", type=Path, help="Path to results_<MxH>_T<T>_cmp_* directory")
    ap.add_argument("--out-prefix", default="", help="Prefix for output filenames")
    args = ap.parse_args()

    result_dir: Path = args.result_dir
    if not result_dir.exists():
        raise SystemExit(f"[ERR] result_dir not found: {result_dir}")

    files = sorted(result_dir.glob("tprime_evolution_*.txt"))
    if not files:
        raise SystemExit(f"[ERR] No tprime_evolution_*.txt under {result_dir}")

    rounds_by_policy = {}
    tighten_by_policy = {}
    for fp in files:
        pol = fp.stem.replace("tprime_evolution_", "")
        _, rounds, tighten = parse_tprime_file(fp)
        rounds_by_policy[pol] = rounds
        tighten_by_policy[pol] = tighten

    # summary
    summary_rows = []
    for pol, d in rounds_by_policy.items():
        rr = [v for v in d["round"] if not math.isnan(v)]
        if not rr:
            continue
        last_idx = len(rr) - 1
        final_Tp = d["T_prime"][last_idx] if last_idx < len(d["T_prime"]) else math.nan
        final_Tt = d["T_target"][last_idx] if last_idx < len(d["T_target"]) else math.nan
        min_maxL = min([v for v in d["maxL"] if not math.isnan(v)]) if d["maxL"] else math.nan
        max_maxL = max([v for v in d["maxL"] if not math.isnan(v)]) if d["maxL"] else math.nan
        meet_rounds = [rr[i] for i, v in enumerate(d["WNS_T"]) if i < len(rr) and not math.isnan(v) and v >= 0]
        first_meet = meet_rounds[0] if meet_rounds else math.nan
        tighten_count = len([r for r in tighten_by_policy.get(pol, {}).get("round", []) if not math.isnan(r)])
        summary_rows.append({
            "policy": pol,
            "rounds": len(rr),
            "tighten_steps": tighten_count,
            "final_T_prime": final_Tp,
            "final_T_target": final_Tt,
            "maxL_min": min_maxL,
            "maxL_max": max_maxL,
            "first_round_meet_T": first_meet,
        })

    out_summary = result_dir / (args.out_prefix + "tprime_summary.csv")
    if summary_rows:
        if HAVE_PANDAS:
            import pandas as pd
            pd.DataFrame(summary_rows).to_csv(out_summary, index=False)
        else:
            with out_summary.open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
                w.writeheader()
                for row in summary_rows:
                    w.writerow(row)

    # plots
    out1 = result_dir / (args.out_prefix + "tprime_vs_round.png")
    out2 = result_dir / (args.out_prefix + "maxL_vs_round.png")
    out3 = result_dir / (args.out_prefix + "wnsT_vs_round.png")

    plot_with_tighten_marks(rounds_by_policy, tighten_by_policy, "T_prime", out1, "T_prime vs. round")
    plot_series_by_policy(rounds_by_policy, "maxL", out2, "maxL vs. round")
    plot_series_by_policy(rounds_by_policy, "WNS_T", out3, "WNS*(T) vs. round")

    print(f"[OK] Wrote summary: {out_summary}")
    print(f"[OK] Wrote figures: {out1}, {out2}, {out3}")


if __name__ == "__main__":
    main()
