"""
Osteotomy planning: a cephalad (proximal) transverse sacral osteotomy placed a fixed
oncologic margin above the most cephalad extent of the tumor, measured perpendicular
to the cut plane.

Plane convention: n . p = d, n is a unit vector pointing cephalad (toward bone that is
KEPT). Everything with n . p < d is the resected specimen (contains the tumor).
"""
import numpy as np


def centroid_curve(bone, dz=2.0, midline_halfwidth=12.0):
    """Centroid of the MIDLINE column (body + canal + crest) of each axial section.
    Whole-slice centroids are biased by the alae, which sit anterior at the S1 level:
    on the TCIA case that tilted the axis by up to 14 deg."""
    zs = np.arange(bone.bounds[0, 2] + 2, bone.bounds[1, 2] - 2, dz)
    pts = []
    for z in zs:
        sec = bone.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
        if sec is None:
            continue
        v = np.asarray(sec.vertices)
        xm = np.median(v[:, 0])
        m = np.abs(v[:, 0] - xm) < midline_halfwidth
        if m.sum() < 6:
            m[:] = True
        pts.append([v[m, 0].mean(), v[m, 1].mean(), z])
    return np.array(pts)


def local_axis(curve, z, window=25.0):  # wide window: robust to 3 mm slice jitter
    sel = np.abs(curve[:, 2] - z) < window
    c = curve[sel]
    # fit y(z) linearly near the level -> sagittal tangent; ignore x drift (midline)
    a = np.polyfit(c[:, 2], c[:, 1], 1)[0]
    t = np.array([0.0, a, 1.0])
    return t / np.linalg.norm(t)


def rot_x(deg):
    a = np.radians(deg)
    return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])


def plan_osteotomy(bone, tumor, margin_mm=10.0, tilt_deg=0.0, iters=4):
    """Return dict with plane normal n, offset d, and diagnostics."""
    curve = centroid_curve(bone)
    tv = np.asarray(tumor.vertices)
    z_level = tv[:, 2].max() + margin_mm
    for _ in range(iters):
        n = rot_x(tilt_deg) @ local_axis(curve, z_level)
        n /= np.linalg.norm(n)
        proj = tv @ n
        d = proj.max() + margin_mm
        # level where plane crosses the bone centroid curve
        f = curve @ n - d
        i = np.argmin(np.abs(f))
        z_level = curve[i, 2]
    q_tumor = tv[np.argmax(proj)]                 # most cephalad tumor point along n
    q_plane = q_tumor + margin_mm * n             # margin-critical point on the plane
    bv = np.asarray(bone.vertices)
    c_rev = curve[::-1]                                   # S1 (top) first
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(c_rev, axis=0), axis=1))])
    level_fraction = float(np.interp(-z_level, -c_rev[:, 2], arc) / arc[-1])
    return {
        "level_fraction": level_fraction,
        "n": n.tolist(), "d": float(d), "margin_mm": margin_mm, "tilt_deg": tilt_deg,
        "q_tumor": q_tumor.tolist(), "q_plane": q_plane.tolist(),
        "z_level_midline": float(z_level),
        "bone_kept_fraction": float(np.mean(bv @ n > d)),
        "min_tumor_to_plane_mm": float(d - proj.max()),
    }
