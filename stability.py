"""
Seating-stability model (linearised rigid-body fit).

Assumption: when the surgeon seats the guide, each contact patch ends up with a random
normal mismatch e_i ~ N(0, tau^2) (residual soft tissue, print error, hand pressure).
The guide settles into the least-squares pose that best explains those mismatches:

    J_i = [ ((p_i - c) x n_i)^T , n_i^T ]        (row of the contact Jacobian)
    xi  = (J^T J)^-1 J^T e        ->   Cov(xi) = tau^2 (J^T J)^-1

xi = (omega, v) is the small rigid twist about the guide centroid c. Contact points are
taken on a `patch_mm` grid so that neighbouring samples are roughly independent
(otherwise sqrt(N) averaging makes every guide look perfect).

Outputs are the cut-plane errors that matter clinically:
  * offset  = displacement along n at the margin-critical point (eats into the margin)
  * angle   = tilt of the cut plane
and the weakest motion mode (the direction the guide can 'slide' most easily).
"""
import numpy as np
import trimesh


def contact_samples(bone, g, patch_mm=5.0):
    step = max(1, int(round(patch_mm / g.params.grid_mm)))
    ii, jj = np.nonzero(g.contact[::step, ::step])
    ii, jj = ii * step, jj * step
    pts = np.column_stack([g.xs[ii], g.surf_y[ii, jj], g.zs[jj]])
    _, _, tri = trimesh.proximity.closest_point(bone, pts)
    nrm = bone.face_normals[tri]
    return pts, nrm


def seating_stability(bone, g, plan, tau=0.5, patch_mm=5.0, n_mc=20000, seed=1):
    n = np.asarray(plan["n"])
    q = np.asarray(plan["q_plane"])
    pts, nrm = contact_samples(bone, g, patch_mm)
    c = pts.mean(0)
    J = np.hstack([np.cross(pts - c, nrm), nrm])
    JTJ = J.T @ J
    # condition number with rotations scaled by a characteristic length
    L = np.sqrt(((pts - c) ** 2).sum(1).mean())
    S = np.diag([1 / L, 1 / L, 1 / L, 1, 1, 1])  # rotations -> mm at radius L
    ev, evec = np.linalg.eigh(S @ JTJ @ S / len(pts))
    cov = tau ** 2 * np.linalg.inv(JTJ)
    rng = np.random.default_rng(seed)
    xi = rng.multivariate_normal(np.zeros(6), cov, n_mc)
    om, v = xi[:, :3], xi[:, 3:]
    disp_q = v + np.cross(om, q - c)
    offset = disp_q @ n
    tilt = np.degrees(np.linalg.norm(om - np.outer(om @ n, n), axis=1))
    weak = evec[:, 0]
    labels = ["rot x (flex/ext)", "rot y (axial)", "rot z (lat. tilt)", "slide x (lateral)",
              "slide y (lift)", "slide z (cranio-caudal)"]
    return {
        "n_contact_patches": int(len(pts)),
        "tau_mm": tau,
        "offset_mean_abs_mm": float(np.mean(np.abs(offset))),
        "offset_p95_mm": float(np.percentile(np.abs(offset), 95)),
        "angle_mean_deg": float(np.mean(tilt)),
        "angle_p95_deg": float(np.percentile(tilt, 95)),
        "eig_min": float(ev[0]), "eig_max": float(ev[-1]),
        "condition": float(ev[-1] / ev[0]),
        "weakest_mode": labels[int(np.argmax(np.abs(weak)))],
        "margin_consumed_p95_pct": float(100 * np.percentile(np.abs(offset), 95) / plan["margin_mm"]),
    }
