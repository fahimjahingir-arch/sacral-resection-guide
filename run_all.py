"""
End-to-end pipeline: CT -> segmentation -> osteotomy plan -> guide variants ->
seating-stability analysis -> FEA -> STL exports + results.json + figures.

  python run_all.py                       # synthetic phantom (pipeline self-test)
  python run_all.py --bone bone.stl --tumor tumor.stl   # your Slicer exports
  python run_all.py --margin 10 --tilt 0 --skip-fea
"""
import argparse
import json
import os
import time

import numpy as np
import trimesh

import guide as G
import plan as PL

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")

VARIANTS = [
    G.GuideParams(name="A_narrow", band_mm=12.0, n_kwires=0),
    G.GuideParams(name="B_wide", band_mm=25.0, n_kwires=0),
    G.GuideParams(name="C_wide_kwires", band_mm=25.0, n_kwires=3),
]


def from_phantom():
    import phantom as P
    import segment as S
    hu, truth = P.make_phantom()
    ix0 = int(round(-P.ORIGIN[0] / P.SPACING))  # x = 0 (midsagittal)
    np.save(os.path.join(OUT, "phantom_hu_midsagittal.npy"), hu[ix0])
    bone_m = S.segment_bone(hu)
    c = np.argwhere(truth["tumor"]).mean(0).round().astype(int)
    tumor_m = S.segment_tumor(hu, c)
    seg = {"dice_bone": S.dice(bone_m, truth["bone"]), "dice_tumor": S.dice(tumor_m, truth["tumor"]),
           "tumor_volume_ml": float(tumor_m.sum() * P.SPACING ** 3 / 1000)}
    bone = S.mask_to_mesh(bone_m, P.SPACING, P.ORIGIN)
    tumor = S.mask_to_mesh(tumor_m, P.SPACING, P.ORIGIN)
    return bone, tumor, seg, hu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bone")
    ap.add_argument("--tumor")
    ap.add_argument("--margin", type=float, default=10.0)
    ap.add_argument("--tilt", type=float, default=0.0)
    ap.add_argument("--tau", type=float, default=0.5, help="seating mismatch RMS, mm")
    ap.add_argument("--skip-fea", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    res = {}
    t0 = time.time()
    if a.bone:
        bone, tumor = trimesh.load(a.bone), trimesh.load(a.tumor)
        res["source"] = {"bone": a.bone, "tumor": a.tumor}
    else:
        bone, tumor, seg, hu = from_phantom()
        res["source"] = "synthetic phantom"
        res["segmentation"] = seg
        print("segmentation", seg)
    bone.export(os.path.join(OUT, "bone.stl"))
    tumor.export(os.path.join(OUT, "tumor.stl"))

    pl = PL.plan_osteotomy(bone, tumor, a.margin, a.tilt)
    res["plan"] = pl
    print("plan: cut at z~%.1f mm, n=%s" % (pl["z_level_midline"], np.round(pl["n"], 3)))

    # resected specimen vs preserved bone (for printing / visualisation)
    n, d = np.asarray(pl["n"]), pl["d"]
    kept_m, spec_m = G.to_manifold(bone).split_by_plane(n.tolist(), d)
    kept, spec = G.from_manifold(kept_m), G.from_manifold(spec_m)
    kept.export(os.path.join(OUT, "bone_preserved.stl"))
    spec.export(os.path.join(OUT, "bone_resected_specimen.stl"))

    import stability as ST
    guides = {}
    res["variants"] = {}
    for p in VARIANTS:
        g = G.build_guide(bone, tumor, pl, p)
        guides[p.name] = g
        g.mesh.export(os.path.join(OUT, f"guide_{p.name}.stl"))
        st = ST.seating_stability(bone, g, pl, tau=a.tau)
        info = {
            "params": G.asdict(p),
            "watertight": bool(g.mesh.is_watertight),
            "volume_cm3": float(g.mesh.volume / 1000),
            "bbox_mm": np.round(g.mesh.extents, 1).tolist(),
            "contact_area_cm2": float(g.contact.sum() * p.grid_mm ** 2 / 100),
            "slot_length_mm": float(g.slot_u_range[1] - g.slot_u_range[0]),
            "kwires": g.kwires,
            "stability": st,
        }
        res["variants"][p.name] = info
        print(p.name, {k: info[k] for k in ("watertight", "volume_cm3", "contact_area_cm2", "slot_length_mm")})
        print("   stability:", {k: st[k] for k in ("offset_p95_mm", "angle_p95_deg", "weakest_mode")})

    if not a.skip_fea:
        import fea
        housing = {}
        for H in (6.0, 10.0, 14.0):
            housing[H] = guides["B_wide"] if H == 10.0 else G.build_guide(
                bone, tumor, pl, G.GuideParams(name=f"B_H{H:g}", band_mm=25.0, housing_height_mm=H))
        res["fea"] = fea.run_study(guides, pl, housing)

    res["runtime_s"] = time.time() - t0
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    import render
    render.all_figures(bone, tumor, guides, pl, res, OUT)
    print("done in %.0f s" % res["runtime_s"])


if __name__ == "__main__":
    main()
