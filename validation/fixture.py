"""
Printable test fixture for TEST 1 / TEST 2.

Two flat datum faces:
  * BACKPLATE (anterior, normal = -v ~ anterior). Lay the fixture on this face while seating the
    guide and sawing. Pressing loads then go straight into the table through a solid anterior
    cradle, with no bending columns. A column-supported design was rejected: 96 mm PLA columns
    would deflect ~1.6 mm under a 20 N press, which is bigger than the effect being measured.
  * BASE (caudal, normal = -n). Stand the fixture on this face to MEASURE: the planned cut
    plane is then horizontal at a known height above the table. With a 1.27 mm shim (or the
    blade) resting in the slot - it sits on the caudal slot wall by gravity - a height gauge at
    4 points gives
        offset = mean(h) - expected_blade_lower_face_mm
        tilt   = plane fit through the 4 points (or a zeroed digital inclinometer)
The cradle has a 5 mm relief groove on the planned plane, so the saw can't cut the fixture,
and both bone pieces stay attached after the cut.

  python fixture.py   -> outputs/test_fixture.stl, outputs/fixture_info.json
"""
import json
import os
import sys

import numpy as np
import trimesh
import manifold3d as m3d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
OUT = os.path.join(os.path.dirname(HERE), "outputs")
import guide as G  # noqa: E402


def main(base_gap=4.0, base_t=8.0, plate_t=12.0, embed=2.0, groove=2.5, grid=1.0):
    res = json.load(open(os.path.join(OUT, "results.json")))
    plan = res["plan"]
    bone = trimesh.load(os.path.join(OUT, "bone.stl"))
    n = np.asarray(plan["n"]); d = plan["d"]
    u, v = G._frame(n)
    R = np.column_stack([u, v, n])            # plan frame (a, b, c) -> world
    local = bone.copy()
    local.vertices = np.asarray(bone.vertices) @ R
    lo, hi = local.bounds

    # anterior surface b_ant(a, c): mirror b so the dorsal-ray helper finds the ANTERIOR face
    mir = local.copy()
    mir.vertices = mir.vertices * np.array([1, -1, 1])
    a_s = np.arange(lo[0] - 2, hi[0] + 2, grid)
    c_s = np.arange(lo[2] - 2, hi[2] + 2, grid)
    negb, _ = G._first_hits(mir, a_s, c_s)
    b_ant = -negb
    b_wall = lo[1] - 3.0
    hit = np.isfinite(b_ant)
    cells = hit[:-1, :-1] & hit[1:, :-1] & hit[:-1, 1:] & hit[1:, 1:]
    cells = G._fix_pinches(cells)
    top = np.where(hit, b_ant + embed, b_wall + 1)
    bottom = np.full_like(top, b_wall - 0.5)
    cradle = G.to_manifold(G.heightmap_solid(a_s, c_s, bottom, top, cells))
    base_top = lo[2] - base_gap
    a0, a1 = lo[0] - 15, hi[0] + 15
    plate = m3d.Manifold.cube([a1 - a0, plate_t, hi[2] - base_top + base_t]).translate(
        [a0, b_wall - plate_t, base_top - base_t])
    base = m3d.Manifold.cube([a1 - a0, hi[1] - b_wall + plate_t + 5, base_t]).translate(
        [a0, b_wall - plate_t, base_top - base_t])
    # cradle reaches down to the base under the apex (caudal specimen stays supported)
    apex_post = m3d.Manifold.cube([16, 14, base_gap + 4]).translate(
        [local.vertices[np.argmin(local.vertices[:, 2]), 0] - 8, b_wall - 1, base_top - 0.5])
    # groove only through the cradle (b > wall), never the backplate, so the halves stay joined
    slab = m3d.Manifold.cube([400, 400, 2 * groove]).translate([-200, b_wall + 1.0, d - groove])
    support = (cradle + plate + base + apex_post) - slab
    fix = support.transform(np.hstack([R, np.zeros((3, 1))]).tolist()) + G.to_manifold(bone)
    parts = fix.decompose()
    # one printable body; tiny negative-volume shells are sealed internal voids -> drop them
    fix = max(parts, key=lambda p: p.volume())
    parts = [p for p in parts if p.volume() > 1.0]
    mesh = G.from_manifold(fix)
    mesh.export(os.path.join(OUT, "test_fixture.stl"))
    gp = res["variants"]["B_wide"]["params"]
    half_w = (gp["blade_mm"] + gp["slot_clearance_mm"]) / 2
    H_plan = d - (base_top - base_t)
    info = {"table_to_planned_plane_mm": float(H_plan),
            "expected_blade_lower_face_mm": float(H_plan - half_w),
            "how_to_use": "seat/cut lying on the BACKPLATE; measure standing on the BASE",
            "note": "the shim rests on the caudal slot wall (gravity), so its lower face sits half a slot width below the plane",
            "n_bodies": len(parts), "watertight": bool(mesh.is_watertight),
            "volume_cm3": float(mesh.volume / 1000), "bbox_mm": np.round(mesh.extents, 1).tolist()}
    json.dump(info, open(os.path.join(OUT, "fixture_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
