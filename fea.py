"""
Linear-elastic FEA of the guide under saw loading (scikit-fem, 8-node hex voxels).

Question answered: how far does the slot wall deflect when the oscillating blade pushes
on it, and what blade-angle error does that produce?

Model
-----
* Geometry: the guide is voxelised directly from its design definition (height-map body,
  housing, slot, bosses, holes) at `pitch` mm.
* Load: F newtons spread over the cephalad slot wall, pushing along -n (into the slot
  wall, i.e. the blade levering against the guide). Linear, so results scale with F.
* Supports:
    'bonded'  - underside contact patches fixed (guide firmly seated / held)
    'pins'    - only the K-wire hole walls fixed (seating contact lost, wires carry it)
    'both'    - contact + pins
* Material: isotropic. Deflection scales ~1/E, so other materials are derived by scaling
  (Poisson effect neglected, <5 %).
Limitations: trilinear hexes are slightly stiff in bending (shear locking) -> deflections
are a lower bound; a pitch convergence check is reported. Bone is rigid. No contact
nonlinearity (a held guide can lift off - 'bonded' is the optimistic case).
"""
import time
import numpy as np
from skfem import MeshHex, Basis, ElementVector, ElementHex1, asm, condense
from skfem.models.elasticity import linear_elasticity, lame_parameters
import pyamg
from scipy.sparse.linalg import cg

import guide as G

MATERIALS = {   # approximate datasheet tensile moduli (MPa) - check your supplier's sheet
    "PA12 (SLS)": 1700.0,
    "Biocompatible SLA resin": 2600.0,
    "Ti-6Al-4V (DMLS)": 110000.0,
}


def inside(g, plan, P):
    p = g.params
    n = np.asarray(plan["n"]); d = plan["d"]
    x, y, z = P.T
    i = np.floor((x - g.xs[0]) / p.grid_mm).astype(int)
    j = np.floor((z - g.zs[0]) / p.grid_mm).astype(int)
    ok = (i >= 0) & (j >= 0) & (i < g.cells.shape[0]) & (j < g.cells.shape[1])
    i, j = np.clip(i, 0, g.cells.shape[0] - 1), np.clip(j, 0, g.cells.shape[1] - 1)
    bot = 0.25 * (g.bottom[i, j] + g.bottom[i + 1, j] + g.bottom[i, j + 1] + g.bottom[i + 1, j + 1])
    top = 0.25 * (g.top[i, j] + g.top[i + 1, j] + g.top[i, j + 1] + g.top[i + 1, j + 1])
    body = ok & g.cells[i, j] & (y > bot) & (y < top)
    for (xk, zk, ytop, yent) in g.kwires:
        r2 = (x - xk) ** 2 + (z - zk) ** 2
        body |= (r2 < p.kwire_boss_r_mm ** 2) & (y > bot) & (y < ytop + p.kwire_boss_h_mm)
        body &= ~(r2 < (p.kwire_hole_mm / 2) ** 2)
    u, v = G._frame(n)
    s = P @ n - d
    us = P @ u
    half_w = (p.blade_mm + p.slot_clearance_mm) / 2
    body &= ~((np.abs(s) < half_w) & (us > g.slot_u_range[0]) & (us < g.slot_u_range[1]))
    return body, bot, top


def build_mesh(g, plan, pitch):
    b = g.mesh.bounds
    ax = [np.arange(b[0, k] - pitch, b[1, k] + 2 * pitch, pitch) for k in range(3)]
    m = MeshHex.init_tensor(*ax)
    cen = m.p[:, m.t].mean(axis=1).T
    keep, _, _ = inside(g, plan, cen)
    m = m.remove_elements(np.flatnonzero(~keep)).remove_unused_nodes()
    return m


def solve(g, plan, pitch=0.7, E=1700.0, nu=0.39, F=20.0, support="bonded"):
    t0 = time.time()
    p = g.params
    n = np.asarray(plan["n"]); d = plan["d"]
    m = build_mesh(g, plan, pitch)
    e = ElementVector(ElementHex1())
    basis = Basis(m, e, intorder=2)
    K = asm(linear_elasticity(*lame_parameters(E, nu)), basis)
    X = m.p.T
    _, bot, top = inside(g, plan, X)
    i = np.clip(np.round((X[:, 0] - g.xs[0]) / p.grid_mm).astype(int), 0, len(g.xs) - 1)
    j = np.clip(np.round((X[:, 2] - g.zs[0]) / p.grid_mm).astype(int), 0, len(g.zs) - 1)
    fixed_nodes = np.zeros(len(X), bool)
    if support in ("bonded", "both"):
        fixed_nodes |= g.contact[i, j] & (X[:, 1] < bot + 0.75 * pitch)
    if support in ("pins", "both"):
        for (xk, zk, ytop, yent) in g.kwires:
            r = np.hypot(X[:, 0] - xk, X[:, 2] - zk)
            fixed_nodes |= r < p.kwire_hole_mm / 2 + 1.01 * pitch
    if fixed_nodes.sum() == 0:
        raise ValueError("no supports for case " + support)
    half_w = (p.blade_mm + p.slot_clearance_mm) / 2
    u_ax, _ = G._frame(n)
    # loaded nodes = mesh-boundary nodes forming the cephalad slot wall. (Selecting by a
    # distance band alone picks interior nodes at some pitches -> mesh-dependent loading.)
    s = X @ n - d
    us = X @ u_ax
    bnd = np.zeros(len(X), bool); bnd[m.boundary_nodes()] = True
    wall = bnd & (s > 0) & (s < half_w + pitch) & (us > g.slot_u_range[0] + pitch) & (us < g.slot_u_range[1] - pitch)
    f = np.zeros(K.shape[0])
    nd = basis.nodal_dofs  # (3, nnodes)
    for k in range(3):
        f[nd[k, wall]] += -F * n[k] / wall.sum()
    D = nd[:, fixed_nodes].ravel()
    Kc, fc, x0, I = condense(K, f, D=D)
    ml = pyamg.smoothed_aggregation_solver(Kc, B=None, max_coarse=500)
    xi, info = cg(Kc, fc, M=ml.aspreconditioner(), rtol=1e-8, maxiter=4000)
    uvec = np.zeros(K.shape[0]); uvec[I] = xi
    U = np.stack([uvec[nd[k]] for k in range(3)], 1)
    un = U[wall] @ n
    lh = X[:, 1] - bot                      # height above local underside
    hmax = lh[wall].max()
    top_edge = wall & (lh > hmax - 1.0)
    root = wall & (lh < 1.0)
    delta_top = float(np.abs(U[top_edge] @ n).mean())
    delta_root = float(np.abs(U[root] @ n).mean()) if root.any() else 0.0
    ang = np.degrees(np.arctan2(abs(delta_top - delta_root), hmax))
    return {
        "support": support, "E_MPa": E, "F_N": F, "pitch_mm": pitch,
        "n_elements": int(m.t.shape[1]), "n_dof": int(K.shape[0]), "cg_info": int(info),
        "wall_max_deflection_mm": float(np.abs(un).max()),
        "wall_mean_deflection_mm": float(np.abs(un).mean()),
        "n_loaded_nodes": int(wall.sum()),
        "wall_top_edge_mm": delta_top, "wall_root_mm": delta_root,
        "wall_height_mm": float(hmax),
        "blade_angle_error_deg": float(ang),
        "max_total_disp_mm": float(np.linalg.norm(U, axis=1).max()),
        "solve_s": round(time.time() - t0, 1),
    }


def run_study(guides, plan, housing_guides=None, F=20.0, pitch=0.8):
    out = {"assumptions": {"F_N": F, "pitch_mm": pitch, "load": "uniform on cephalad slot wall, along -n",
                           "materials_MPa": MATERIALS}}
    E0, nu0 = MATERIALS["PA12 (SLS)"], 0.39
    cases = []
    for name, g in guides.items():
        sups = ["bonded"] + (["pins", "both"] if g.kwires else [])
        for s in sups:
            r = solve(g, plan, pitch, E0, nu0, F, s)
            r["variant"] = name
            cases.append(r)
            print("FEA", name, s, {k: round(r[k], 4) for k in ("wall_max_deflection_mm", "blade_angle_error_deg")},
                  r["n_elements"], "el", r["solve_s"], "s")
    out["variants_PA12"] = cases
    # material scaling on variant B bonded
    base = next(c for c in cases if c["variant"] == "B_wide" and c["support"] == "bonded")
    out["materials_B_bonded"] = {
        mat: {"wall_max_deflection_mm": base["wall_max_deflection_mm"] * E0 / E,
              "blade_angle_error_deg": base["blade_angle_error_deg"] * E0 / E}
        for mat, E in MATERIALS.items()}
    if housing_guides:
        hs = []
        for H, g in housing_guides.items():
            r = solve(g, plan, pitch, E0, nu0, F, "bonded")
            r["housing_height_mm"] = H
            hs.append(r)
            print("FEA housing", H, round(r["wall_max_deflection_mm"], 4), round(r["blade_angle_error_deg"], 4))
        out["housing_height_B_PA12"] = hs
    # mesh convergence on B bonded
    conv = []
    for pt in (1.2, 1.0, 0.8, 0.6):
        r = base if pt == pitch else solve(guides["B_wide"], plan, pt, E0, nu0, F, "bonded")
        conv.append({"pitch_mm": pt, "n_elements": r["n_elements"],
                     "wall_max_deflection_mm": r["wall_max_deflection_mm"]})
        print("FEA convergence", conv[-1])
    out["convergence_B"] = conv
    out["slot_play"] = slot_play(guides["B_wide"], housing_guides)
    return out


def slot_play(g, housing_guides=None):
    """Rigid-body blade wobble allowed by slot clearance: atan(clearance / guided depth)."""
    p = g.params
    rows = []
    gs = housing_guides or {p.housing_height_mm: g}
    for H, gg in gs.items():
        depth = gg.params.base_thickness_mm + H
        rows.append({"housing_height_mm": H, "guided_depth_mm": depth,
                     "clearance_mm": p.slot_clearance_mm,
                     "play_angle_deg": float(np.degrees(np.arctan(p.slot_clearance_mm / depth)))})
    return rows
