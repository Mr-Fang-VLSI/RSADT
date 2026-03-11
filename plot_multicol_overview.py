#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize multi-column (K=2/4) results produced by run_multicol_K2K4.sh.
Generates:
  - Global overlay plots (K=2 vs K=4): runtime, HPWL_eq, HPWL_w/HPWL_eq, first_meet_round, tighten_steps, avg_solve_ms vs T
  - Per-K overlays across T: T_prime vs round, maxL vs round (by parsing logs)
  - Per-K overview CSVs; optional Excel aggregation if pandas+openpyxl available.
Usage:
  python3 plot_multicol_overview.py <RESULT_ROOT> [--out-prefix PFX]
"""
import argparse, csv, math, re
from pathlib import Path
from typing import Dict, List

try:
    import pandas as pd
    HAVE_PANDAS = True
except Exception:
    HAVE_PANDAS = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_csv_rows(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({(k or '').strip(): (v or '').strip() for k,v in r.items()})
    return rows

def to_num(s: str):
    if s is None: return math.nan
    t = s.strip()
    if t=='' or t.lower()=='na': return math.nan
    try:
        if any(c in t for c in ('.','e','E')): return float(t)
        return int(t)
    except Exception:
        try: return float(t)
        except Exception: return math.nan

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def parse_round_series_from_log(fp: Path):
    d = {"round": [], "T_prime": [], "maxL": [], "HPWL_eq": [], "HPWL_w": [], "WNS_T": [], "solve_ms": []}
    mT = re.search(r"_T_([0-9]+)\.log$", fp.name)
    T_val = int(mT.group(1)) if mT else None
    with fp.open('r', encoding='utf-8', errors='ignore') as f:
        for raw in f:
            line = raw.strip()
            if '[Round' not in line: continue
            m = re.search(r"\[Round\s+([0-9]+)\]", line); r = int(m.group(1)) if m else None
            m = re.search(r"T_prime=([0-9]+)", line); Tp = float(m.group(1)) if m else None
            m = re.search(r"T_target=([0-9]+)", line); Tt = float(m.group(1)) if m else None
            m = re.search(r"maxL=([0-9]+)", line); ml = float(m.group(1)) if m else None
            m = re.search(r"HPWL_eq=([0-9]+)", line); he = float(m.group(1)) if m else None
            m = re.search(r"HPWL_w=([0-9]+)", line); hw = float(m.group(1)) if m else None
            m = re.search(r"solve_ms=([0-9]*\.?[0-9]+)", line); ms = float(m.group(1)) if m else None
            m = re.search(r"WNS\*\(T\)=(-?[0-9]+)", line); wT = float(m.group(1)) if m else None
            if r is None: continue
            d['round'].append(float(r))
            d['T_prime'].append(Tp if Tp is not None else math.nan)
            d['maxL'].append(ml if ml is not None else math.nan)
            d['HPWL_eq'].append(he if he is not None else math.nan)
            d['HPWL_w'].append(hw if hw is not None else math.nan)
            d['solve_ms'].append(ms if ms is not None else math.nan)
            if wT is not None:
                d['WNS_T'].append(wT)
            else:
                T_use = Tt if Tt is not None else (float(T_val) if T_val is not None else math.nan)
                d['WNS_T'].append(T_use - ml if (not math.isnan(T_use) and ml is not None) else math.nan)
    return d

def plot_overlay_byK(byK, xkey, ykey, title, out_png):
    plt.figure()
    for K in sorted(byK.keys()):
        x = byK[K].get(xkey, []); y = byK[K].get(ykey, [])
        x2=[]; y2=[]
        for i, xv in enumerate(x):
            if math.isnan(xv): continue
            yv = y[i] if i < len(y) else math.nan
            if math.isnan(yv): continue
            x2.append(xv); y2.append(yv)
        if x2 and y2: plt.plot(x2, y2, label=f"K={K}")
    plt.xlabel(xkey); plt.ylabel(ykey); plt.title(title); plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5); plt.tight_layout(); plt.savefig(out_png); plt.close()

def process_root(root: Path, out_prefix: str=''):
    glob = root / 'summary_multicol_all.txt'
    if not glob.exists(): raise SystemExit(f"[ERR] not found: {glob}")
    rows = read_csv_rows(glob)

    # group by K
    byK = {}
    for r in rows:
        try: K = int(r.get('K','0') or '0')
        except: continue
        if K not in (2,4): continue
        d = byK.setdefault(K, {'T':[], 'total_runtime_sec':[], 'last_HPWL_eq':[],
                               'last_HPWL_w':[], 'first_meet_T_round':[],
                               'tighten_steps':[], 'avg_solve_ms':[]})
        d['T'].append(to_num(r.get('T','')))
        d['total_runtime_sec'].append(to_num(r.get('total_runtime_sec','')))
        d['last_HPWL_eq'].append(to_num(r.get('last_HPWL_eq','')))
        d['last_HPWL_w'].append(to_num(r.get('last_HPWL_w','')))
        d['first_meet_T_round'].append(to_num(r.get('first_meet_T_round','')))
        d['tighten_steps'].append(to_num(r.get('tighten_steps','')))
        d['avg_solve_ms'].append(to_num(r.get('avg_solve_ms','')))

    # sort by T and compute ratio
    for K, d in list(byK.items()):
        zipped = list(zip(d['T'], d['total_runtime_sec'], d['last_HPWL_eq'], d['last_HPWL_w'],
                          d['first_meet_T_round'], d['tighten_steps'], d['avg_solve_ms']))
        zipped.sort(key=lambda t: (math.inf if math.isnan(t[0]) else t[0]))
        T_sorted  = [z[0] for z in zipped]
        runtime   = [z[1] for z in zipped]
        hpwl_eq   = [z[2] for z in zipped]
        hpwl_w    = [z[3] for z in zipped]
        firstmeet = [z[4] for z in zipped]
        tighten   = [z[5] for z in zipped]
        avgms     = [z[6] for z in zipped]
        ratio = []
        for he, hw in zip(hpwl_eq, hpwl_w):
            if he and (not math.isnan(he)) and he != 0 and hw and (not math.isnan(hw)):
                ratio.append(hw/he)
            else:
                ratio.append(math.nan)
        byK[K] = {'T':T_sorted, 'total_runtime_sec':runtime, 'last_HPWL_eq':hpwl_eq,
                  'HPWL_ratio':ratio, 'first_meet_T_round':firstmeet,
                  'tighten_steps':tighten, 'avg_solve_ms':avgms}

    plots_dir = root / 'plots'; ensure_dir(plots_dir)
    plot_overlay_byK(byK, 'T', 'total_runtime_sec', 'Runtime vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_T_vs_runtime_byK.png')
    plot_overlay_byK(byK, 'T', 'last_HPWL_eq', 'HPWL_eq vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_T_vs_HPWL_eq_byK.png')
    plot_overlay_byK(byK, 'T', 'HPWL_ratio', 'HPWL_w / HPWL_eq vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_HPWL_ratio_vs_T_byK.png')
    plot_overlay_byK(byK, 'T', 'first_meet_T_round', 'First meet T (round) vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_first_meet_round_vs_T_byK.png')
    plot_overlay_byK(byK, 'T', 'tighten_steps', 'Tighten steps vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_tighten_steps_vs_T_byK.png')
    plot_overlay_byK(byK, 'T', 'avg_solve_ms', 'Avg DP solve_ms vs T (K=2 vs K=4)', plots_dir / f'{out_prefix}global_avg_solve_ms_vs_T_byK.png')

    # Per-K round overlays
    for kd in [p for p in root.iterdir() if p.is_dir() and re.match(r'K[24]_\d+x\d+_T', p.name)]:
        logs = sorted((kd/'logs').glob('*.log'))
        series_by_T = {}
        for lg in logs:
            m = re.search(r'_T_([0-9]+)\.log$', lg.name)
            if not m: continue
            Tval = int(m.group(1))
            series_by_T[Tval] = parse_round_series_from_log(lg)
        if not series_by_T: continue

        kplots = kd / 'plots'; ensure_dir(kplots)
        # T_prime vs round
        plt.figure()
        for Tval in sorted(series_by_T.keys()):
            d = series_by_T[Tval]
            rr = [v for v in d['round'] if not math.isnan(v)]
            yy = [d['T_prime'][i] for i,_ in enumerate(d['round']) if not math.isnan(d['round'][i])]
            if rr and yy: plt.plot(rr, yy, label=f'T={Tval}')
        plt.xlabel('round'); plt.ylabel('T_prime'); plt.title(f'T_prime vs round ({kd.name})')
        plt.legend(); plt.grid(True, linestyle='--', linewidth=0.5); plt.tight_layout()
        plt.savefig(kplots / f'{out_prefix}Tprime_vs_round_allT_{kd.name}.png'); plt.close()

        # maxL vs round
        plt.figure()
        for Tval in sorted(series_by_T.keys()):
            d = series_by_T[Tval]
            rr = [v for v in d['round'] if not math.isnan(v)]
            yy = [d['maxL'][i] for i,_ in enumerate(d['round']) if not math.isnan(d['round'][i])]
            if rr and yy: plt.plot(rr, yy, label=f'T={Tval}')
        plt.xlabel('round'); plt.ylabel('maxL'); plt.title(f'maxL vs round ({kd.name})')
        plt.legend(); plt.grid(True, linestyle='--', linewidth=0.5); plt.tight_layout()
        plt.savefig(kplots / f'{out_prefix}maxL_vs_round_allT_{kd.name}.png'); plt.close()

        # per-K short overview CSV
        ov = kd / f'{out_prefix}overview_{kd.name}.csv'
        with ov.open('w', newline='', encoding='utf-8') as f:
            w = csv.writer(f); w.writerow(['T','round_len','final_T_prime','final_maxL','final_HPWL_eq','final_HPWL_w'])
            for Tval in sorted(series_by_T.keys()):
                d = series_by_T[Tval]; rr = [v for v in d['round'] if not math.isnan(v)]
                idx = len(rr)-1 if rr else -1
                fTp = d['T_prime'][idx] if idx>=0 else ''
                fL  = d['maxL'][idx] if idx>=0 else ''
                fHe = d['HPWL_eq'][idx] if idx>=0 else ''
                fHw = d['HPWL_w'][idx] if idx>=0 else ''
                w.writerow([Tval, len(rr), fTp, fL, fHe, fHw])

    # Excel (optional)
    if HAVE_PANDAS:
        try:
            xlsx = root / (out_prefix + 'multicol_overview.xlsx')
            with pd.ExcelWriter(xlsx, engine='openpyxl') as writer:
                for K, d in byK.items():
                    import pandas as pd
                    df = pd.DataFrame({
                        'T': d['T'],
                        'runtime_sec': d['total_runtime_sec'],
                        'HPWL_eq': d['last_HPWL_eq'],
                        'HPWL_ratio': d['HPWL_ratio'],
                        'first_meet_round': d['first_meet_T_round'],
                        'tighten_steps': d['tighten_steps'],
                        'avg_solve_ms': d['avg_solve_ms'],
                    })
                    df.to_excel(writer, sheet_name=f'global_K{K}', index=False)
                for kd in [p for p in root.iterdir() if p.is_dir() and re.match(r'K[24]_\d+x\d+_T', p.name)]:
                    ov = kd / f'{out_prefix}overview_{kd.name}.csv'
                    if ov.exists():
                        df = pd.read_csv(ov); df.to_excel(writer, sheet_name=f'{kd.name}', index=False)
            print(f'[OK] Wrote Excel: {xlsx}')
        except Exception as e:
            print(f'[WARN] Excel export skipped: {e}')

    print(f'[OK] Plots saved under: {root}/plots')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('result_root', type=Path, help='Path to results_multicol_<MxH>_<timestamp>')
    ap.add_argument('--out-prefix', default='', help='Prefix for output filenames')
    args = ap.parse_args()
    if not args.result_root.exists():
        raise SystemExit(f'[ERR] not found: {args.result_root}')
    process_root(args.result_root, args.out_prefix)

if __name__ == '__main__':
    main()
