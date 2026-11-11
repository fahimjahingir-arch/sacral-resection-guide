"""Re-draw all figures from outputs/ without re-running FEA (guides are rebuilt, ~10 s)."""
import json
import os
import sys

import trimesh

import guide as G
import render
from run_all import VARIANTS, OUT

out = sys.argv[1] if len(sys.argv) > 1 else OUT
res = json.load(open(os.path.join(OUT, "results.json")))
bone = trimesh.load(os.path.join(OUT, "bone.stl"))
tumor = trimesh.load(os.path.join(OUT, "tumor.stl"))
guides = {p.name: G.build_guide(bone, tumor, res["plan"], p) for p in VARIANTS}
os.makedirs(out, exist_ok=True)
if out != OUT and os.path.exists(os.path.join(OUT, "phantom_hu_midsagittal.npy")):
    import shutil
    shutil.copy(os.path.join(OUT, "phantom_hu_midsagittal.npy"), out)
render.all_figures(bone, tumor, guides, res["plan"], res, out)
print("figures ->", os.path.join(out, "figures"))
