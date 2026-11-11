"""
Sterilisation dimensional-stability test (TEST 3 in PROTOCOLS.md).

  python autoclave_coupon.py --make      # writes coupon STL + empty data sheet
  python autoclave_coupon.py data.csv    # analyse pre/post measurements

The coupon reproduces the guide's critical features at the same sizes:
  * saw slot 1.40 mm wide x 20 mm long x 14 mm deep (housing)  -> blade fit
  * 3 x 1.8 mm K-wire holes                                     -> wire fit
  * 60 mm long thin rib (3 mm)                                  -> warpage / bow
  * 20.00 mm reference block length                             -> bulk shrinkage

CSV columns: material, coupon, feature, nominal_mm, pre_mm, post_mm, cycles
"""
import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
OUT = os.path.join(os.path.dirname(HERE), "outputs")
BLADE = 1.27
FEATURES = [("slot_width", 1.40), ("hole1_dia", 1.80), ("hole2_dia", 1.80), ("hole3_dia", 1.80),
            ("block_length", 20.00), ("rib_length", 60.00), ("rib_bow", 0.00), ("slot_depth", 14.00)]


def make():
    import manifold3d as m3d
    import guide as G
    C = m3d.Manifold.cube
    body = C([20, 30, 14])                                   # block_length = 20 mm in x
    slot = C([1.40, 20.0, 15]).translate([9.30, 5.0, -0.5])  # through-slot, closed ends
    rib = C([60, 3, 6]).translate([20, 25, 0])               # cantilever rib for bow
    plate = C([20, 20, 4]).translate([20, 0, 0])
    holes = m3d.Manifold()
    for k in range(3):
        holes = holes + m3d.Manifold.cylinder(6, 0.9, 0.9, 48).translate([24 + 6 * k, 10, -1])
    part = (body - slot) + rib + (plate - holes)
    mesh = G.from_manifold(part)
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, "autoclave_coupon.stl")
    mesh.export(p)
    d = os.path.join(HERE, "autoclave_datasheet.csv")
    with open(d, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["material", "coupon", "feature", "nominal_mm", "pre_mm", "post_mm", "cycles"])
        for mat in ("PA12_SLS", "SLA_bio_resin"):
            for c in range(1, 6):
                for feat, nom in FEATURES:
                    w.writerow([mat, c, feat, nom, "", "", 1])
    print("wrote", p, "and", d, "| volume %.1f cm3" % (mesh.volume / 1000))


def analyze(path):
    rows = [r for r in csv.DictReader(open(path, newline="")) if r["pre_mm"] and r["post_mm"]]
    mats = sorted({r["material"] for r in rows})
    for m in mats:
        print(f"\n== {m} ==")
        for feat, nom in FEATURES:
            R = [r for r in rows if r["material"] == m and r["feature"] == feat]
            if not R:
                continue
            pre = np.array([float(r["pre_mm"]) for r in R]); post = np.array([float(r["post_mm"]) for r in R])
            dlt = post - pre
            pct = 100 * dlt / nom if nom else np.full_like(dlt, np.nan)
            line = (f"{feat:13s} n={len(R)}  pre {pre.mean():.3f}  post {post.mean():.3f}  "
                    f"Δ {dlt.mean():+.3f}±{dlt.std(ddof=1) if len(R) > 1 else 0:.3f} mm")
            if nom:
                line += f"  ({np.nanmean(pct):+.2f} %)"
            if feat == "slot_width":
                fit = post - BLADE
                line += f"  | blade clearance after: {fit.min():.3f}..{fit.max():.3f} mm" + \
                        ("  ** BLADE MAY BIND **" if fit.min() < 0.05 else "")
            print(line)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--make", action="store_true")
    a = ap.parse_args()
    make() if a.make else analyze(a.csv)
