"""
Cut-accuracy analysis for the physical fit test (TEST 1 in PROTOCOLS.md).

Input : a 3D scan (STL) of the PRESERVED (cephalad) piece of a printed sacrum after you
        cut it through the guide, plus the planned bone model and plan (results.json).
Steps : 1) ICP-register the scan to the planned bone using surfaces away from the cut
        2) pick the scanned cut face (normals ~ plan normal, near the plane)
        3) least-squares plane fit -> compare with the planned plane
Output: offset at the margin-critical point (mm, + = cut moved cephalad, i.e. extra
        bone removed; - = cut moved TOWARD the tumor, i.e. margin lost), angle (deg),
        cut-face flatness (RMS, mm).

  python analyze_cut.py --scan specimen01.stl
  python analyze_cut.py --selftest        # synthetic check with a known error
"""
import argparse
import json
import os
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
OUT = os.path.join(os.path.dirname(HERE), "outputs")


def fit_plane(P):
    c = P.mean(0)
    _, s, vt = np.linalg.svd(P - c)
    nrm = vt[2]
    rms = s[2] / np.sqrt(len(P))
    return c, nrm, rms


def icp_trimmed(src, ref_pts, ref_nrm, init, iters=60):
    """Point-to-plane ICP with distance trimming. Scan points that have no counterpart on
    the intact model (the new cut face, scan edges) are rejected instead of biasing the fit."""
    from scipy.spatial import cKDTree
    tree = cKDTree(ref_pts)
    T = init.copy()
    for k in range(iters):
        P = src @ T[:3, :3].T + T[:3, 3]
        dist, idx = tree.query(P)
        thr = 5.0 if k < 10 else max(0.4, 3.0 * np.median(dist))
        keep = dist < thr
        p, q, nn = P[keep], ref_pts[idx[keep]], ref_nrm[idx[keep]]
        r = np.einsum("ij,ij->i", p - q, nn)
        A = np.hstack([np.cross(p, nn), nn])
        x, *_ = np.linalg.lstsq(A, -r, rcond=None)
        w, t = x[:3], x[3:]
        th = np.linalg.norm(w)
        Rk = trimesh.transformations.rotation_matrix(th, w / th)[:3, :3] if th > 1e-12 else np.eye(3)
        D = np.eye(4); D[:3, :3] = Rk; D[:3, 3] = t
        T = D @ T
        if th < 1e-7 and np.linalg.norm(t) < 1e-6:
            break
    rms = float(np.sqrt(np.mean(r ** 2)))
    return T, rms, float(keep.mean())


def analyze(scan, bone, plan, init=None):
    n = np.asarray(plan["n"]); d = plan["d"]; q = np.asarray(plan["q_plane"])
    sp = scan.sample(20000)
    # 1) register the scan to the intact planned model (dense samples with normals)
    ref, fi = trimesh.sample.sample_surface(bone, 300000)
    T, cost, inlier = icp_trimmed(sp, ref, bone.face_normals[fi], np.eye(4) if init is None else init)
    reg = scan.copy()
    reg.apply_transform(T)
    # 2) cut-face triangles: normal anti-parallel to n (the preserved piece's cut face looks caudal)
    fn = reg.face_normals
    cen = reg.triangles_center
    sel = (fn @ -n > 0.9) & (np.abs(cen @ n - d) < 6.0)
    pts, fidx = trimesh.sample.sample_surface(reg.submesh([np.flatnonzero(sel)], append=True), 8000)
    c, m, rms = fit_plane(pts)
    if m @ n < 0:
        m = -m
    dm = m @ c
    # 3) compare: offset of the fitted plane at the margin-critical point, along n
    t = (dm - m @ q) / (m @ n)           # q + t n lies on the fitted plane
    angle = np.degrees(np.arccos(np.clip(abs(m @ n), -1, 1)))
    return {"offset_mm": float(t), "angle_deg": float(angle), "flatness_rms_mm": float(rms),
            "icp_rms_mm": float(cost), "icp_inlier_fraction": inlier, "n_cut_points": int(len(pts)),
            "margin_remaining_mm": float(plan["margin_mm"] + t)}


def _rot(axis, deg):
    return trimesh.transformations.rotation_matrix(np.radians(deg), axis)


def selftest():
    import guide as G
    res = json.load(open(os.path.join(OUT, "results.json")))
    plan = res["plan"]
    bone = trimesh.load(os.path.join(OUT, "bone.stl"))
    n = np.asarray(plan["n"]); d = plan["d"]; q = np.asarray(plan["q_plane"])
    rows = []
    rng = np.random.default_rng(3)
    for true_off, true_ang in [(0.0, 0.0), (-1.5, 0.0), (1.0, 3.0), (-0.8, 5.0)]:
        # simulated actual cut plane: rotate n about the lateral axis through q, then shift
        R = _rot([1, 0, 0], true_ang)[:3, :3]
        n2 = R @ n
        d2 = n2 @ (q + true_off * n)
        kept, _ = G.to_manifold(bone).split_by_plane(n2.tolist(), float(d2))
        piece = G.from_manifold(kept)
        piece.vertices += rng.normal(0, 0.08, piece.vertices.shape)    # scanner noise
        mis = _rot([0.3, 1, 0.2], 4.0)
        mis[:3, 3] = [3.0, -2.0, 4.0]                                  # pose in scanner frame
        piece.apply_transform(mis)
        r = analyze(piece, bone, plan, init=np.linalg.inv(mis))
        # true offset measured at q along n (same definition as analyze)
        t_true = (d2 - n2 @ q) / (n2 @ n)
        rows.append({"true_offset_mm": round(float(t_true), 3), "true_angle_deg": true_ang,
                     **{k: round(v, 3) for k, v in r.items()}})
        print(rows[-1])
    err_o = max(abs(r["true_offset_mm"] - r["offset_mm"]) for r in rows)
    err_a = max(abs(r["true_angle_deg"] - r["angle_deg"]) for r in rows)
    print("self-test max error: offset %.3f mm, angle %.3f deg" % (err_o, err_a))
    json.dump({"cases": rows, "max_offset_error_mm": err_o, "max_angle_error_deg": err_a},
              open(os.path.join(OUT, "selftest_analyze_cut.json"), "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan")
    ap.add_argument("--results", default=os.path.join(OUT, "results.json"))
    ap.add_argument("--bone", default=os.path.join(OUT, "bone.stl"))
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        plan = json.load(open(a.results))["plan"]
        print(json.dumps(analyze(trimesh.load(a.scan), trimesh.load(a.bone), plan), indent=1))
