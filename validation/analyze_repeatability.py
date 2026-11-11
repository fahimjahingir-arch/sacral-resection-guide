"""
Placement repeatability (TEST 2 in PROTOCOLS.md).

CSV columns: variant, operator, trial, offset_mm, angle_deg
(offset/angle per placement, measured with the slot-gauge method or analyze_cut.py)

Reports per variant: mean, SD, 95 % limits, within-operator (repeatability) SD,
between-operator (reproducibility) SD from a one-way random-effects ANOVA, and a
Welch t-test / Levene test comparing variances between variants.

  python analyze_repeatability.py data.csv
  python analyze_repeatability.py --make-template   # writes an empty data sheet
  python analyze_repeatability.py --example         # SYNTHETIC numbers, demo only
"""
import argparse
import csv
import itertools
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def variance_components(groups):
    k = len(groups)
    n = np.array([len(g) for g in groups])
    allv = np.concatenate(groups)
    grand = allv.mean()
    ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
    ssb = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    msw = ssw / (n.sum() - k)
    msb = ssb / (k - 1) if k > 1 else 0.0
    n0 = (n.sum() - (n ** 2).sum() / n.sum()) / (k - 1) if k > 1 else n[0]
    var_b = max(0.0, (msb - msw) / n0)
    return np.sqrt(msw), np.sqrt(var_b)


def analyze(rows):
    out = {}
    variants = sorted({r["variant"] for r in rows})
    for v in variants:
        R = [r for r in rows if r["variant"] == v]
        res = {"n": len(R), "operators": len({r["operator"] for r in R})}
        for key in ("offset_mm", "angle_deg"):
            x = np.array([float(r[key]) for r in R])
            ops = sorted({r["operator"] for r in R})
            groups = [np.array([float(r[key]) for r in R if r["operator"] == o]) for o in ops]
            sw, sb = variance_components(groups)
            res[key] = {"mean": x.mean(), "sd": x.std(ddof=1), "lo95": x.mean() - 1.96 * x.std(ddof=1),
                        "hi95": x.mean() + 1.96 * x.std(ddof=1), "abs_p95": np.percentile(np.abs(x), 95),
                        "within_operator_sd": sw, "between_operator_sd": sb}
        out[v] = res
    comps = {}
    for a, b in itertools.combinations(variants, 2):
        xa = np.array([float(r["offset_mm"]) for r in rows if r["variant"] == a])
        xb = np.array([float(r["offset_mm"]) for r in rows if r["variant"] == b])
        comps[f"{a} vs {b}"] = {"levene_p_offset_variance": float(stats.levene(xa, xb).pvalue),
                                "welch_p_offset_mean": float(stats.ttest_ind(xa, xb, equal_var=False).pvalue)}
    return out, comps


def load(path):
    with open(path, newline="") as f:
        return [r for r in csv.DictReader(f) if r.get("offset_mm")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--make-template", action="store_true")
    ap.add_argument("--example", action="store_true")
    a = ap.parse_args()
    if a.make_template:
        p = os.path.join(HERE, "repeatability_datasheet.csv")
        with open(p, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["variant", "operator", "trial", "offset_mm", "angle_deg", "notes"])
            for v in ("A_narrow", "B_wide", "C_wide_kwires"):
                for o in ("op1", "op2", "op3", "op4", "op5"):
                    for t in range(1, 11):
                        w.writerow([v, o, t, "", "", ""])
        print("wrote", p)
        return
    if a.example:
        rng = np.random.default_rng(7)
        rows = []
        for v, sd in (("A_narrow", 0.9), ("B_wide", 0.5), ("C_wide_kwires", 0.45)):
            for o in range(5):
                bias = rng.normal(0, sd / 2)
                for t in range(10):
                    rows.append({"variant": v, "operator": f"op{o+1}", "trial": t + 1,
                                 "offset_mm": bias + rng.normal(0, sd), "angle_deg": abs(rng.normal(0, 2 * sd))})
        print("*** SYNTHETIC EXAMPLE DATA - NOT A MEASUREMENT ***")
    else:
        if not a.csv:
            sys.exit("give a CSV, --make-template or --example")
        rows = load(a.csv)
    res, comps = analyze(rows)
    for v, r in res.items():
        o, g = r["offset_mm"], r["angle_deg"]
        print(f"{v:15s} n={r['n']:3d}  offset {o['mean']:+.2f}±{o['sd']:.2f} mm "
              f"(95% LoA {o['lo95']:+.2f}..{o['hi95']:+.2f}; within-op SD {o['within_operator_sd']:.2f}, "
              f"between-op SD {o['between_operator_sd']:.2f})  angle {g['mean']:.2f}±{g['sd']:.2f}°")
    for k, c in comps.items():
        print(f"  {k}: Levene p={c['levene_p_offset_variance']:.3g}, Welch p={c['welch_p_offset_mean']:.3g}")


if __name__ == "__main__":
    main()
