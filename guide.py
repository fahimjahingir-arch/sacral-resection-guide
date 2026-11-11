"""
Surface-matched cutting-guide generator.

Design rules encoded here
-------------------------
* Insertion direction = -y (posterior approach, guide lowered straight onto the dorsal
  sacrum). The guide underside is a height-map of the bone seen from that direction, so
  by construction it has NO undercuts and seats along a single path.
* Underside = bone surface + clearance (accounts for print tolerance + residual
  periosteum). Foramina / notches are bridged (no contact, +relief) instead of dipping in.
* No contact within `tumor_keepout_mm` of tumor (avoid contamination / seeding).
* Slot = blade thickness + clearance, aligned EXACTLY to the planned plane, with a raised
  housing for blade guidance depth and closed end-bridges so the guide stays one piece.
* Optional K-wire holes parallel to the insertion direction (guide can be lifted off over
  the wires), only where the bone is thick enough that the wire cannot breach the
  anterior cortex or enter the canal/foramina, and away from the tumor.
"""
from dataclasses import dataclass, field, asdict
import numpy as np
import trimesh
import manifold3d as m3d
from scipy import ndimage as ndi
from scipy.spatial import cKDTree


@dataclass
class GuideParams:
    name: str = "B_wide"
    band_mm: float = 25.0          # footprint half-height either side of the cut line
    lateral_mm: float = 60.0       # max |x| of footprint
    grid_mm: float = 0.5
    clearance_mm: float = 0.3
    hole_relief_mm: float = 1.0
    base_thickness_mm: float = 4.0
    housing_half_mm: float = 6.0   # housing extends this far each side of the slot
    housing_height_mm: float = 10.0
    blade_mm: float = 1.27         # common sagittal-saw blade thickness
    slot_clearance_mm: float = 0.13
    end_bridge_mm: float = 5.0
    tumor_keepout_mm: float = 3.0
    n_kwires: int = 0
    kwire_hole_mm: float = 1.8     # for 1.6 mm K-wires
    kwire_boss_r_mm: float = 3.5
    kwire_boss_h_mm: float = 8.0
    kwire_depth_mm: float = 15.0
    kwire_safety_mm: float = 3.0


@dataclass
class Guide:
    params: GuideParams
    xs: np.ndarray
    zs: np.ndarray
    bottom: np.ndarray            # (nx, nz) y of underside at grid vertices
    top: np.ndarray
    cells: np.ndarray             # (nx-1, nz-1) bool cell mask
    contact: np.ndarray           # (nx, nz) bool vertex in real bone contact
    surf_y: np.ndarray            # bone surface y (NaN where no bone)
    h: np.ndarray                 # signed distance of bone surface point to plane
    slot_u_range: tuple = None
    kwires: list = field(default_factory=list)  # [(x, z, y_top, y_entry)]
    mesh: trimesh.Trimesh = None


def _first_hits(mesh, xs, zs, y0=300.0):
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    orig = np.column_stack([X.ravel(), np.full(X.size, y0), Z.ravel()])
    dirs = np.tile([0.0, -1.0, 0.0], (X.size, 1))
    loc, ray, _ = mesh.ray.intersects_location(orig, dirs, multiple_hits=True)
    hits = [[] for _ in range(X.size)]
    for p, r in zip(loc, ray):
        hits[r].append(p[1])
    ys = np.full(X.size, np.nan)
    for i, hl in enumerate(hits):
        if hl:
            ys[i] = max(hl)
    return ys.reshape(X.shape), [sorted(h, reverse=True) for h in hits]


def _fill_nan(a, sigma=4.0):
    m = np.isfinite(a)
    num = ndi.gaussian_filter(np.where(m, a, 0.0), sigma)
    den = ndi.gaussian_filter(m.astype(float), sigma)
    out = a.copy()
    out[~m] = (num / np.maximum(den, 1e-9))[~m]
    return out


def _fix_pinches(cells):
    c = cells.copy()
    for _ in range(50):
        a, b = c[:-1, :-1], c[1:, :-1]
        cc, d = c[:-1, 1:], c[1:, 1:]
        p1 = a & d & ~b & ~cc
        p2 = b & cc & ~a & ~d
        if not (p1.any() or p2.any()):
            break
        ii, jj = np.nonzero(p1)
        c[ii, jj] = False
        ii, jj = np.nonzero(p2)
        c[ii + 1, jj] = False
    return c


def heightmap_solid(xs, zs, bottom, top, cells):
    nx, nz = len(xs), len(zs)
    vid = -np.ones((nx, nz), int)
    used = np.zeros((nx, nz), bool)
    ci, cj = np.nonzero(cells)
    for di in (0, 1):
        for dj in (0, 1):
            used[ci + di, cj + dj] = True
    ui, uj = np.nonzero(used)
    vid[ui, uj] = np.arange(len(ui))
    nv = len(ui)
    X, Z = xs[ui], zs[uj]
    V = np.vstack([np.column_stack([X, top[ui, uj], Z]),
                   np.column_stack([X, bottom[ui, uj], Z])])
    v00, v10 = vid[ci, cj], vid[ci + 1, cj]
    v01, v11 = vid[ci, cj + 1], vid[ci + 1, cj + 1]
    F = [np.column_stack([v00, v11, v10]), np.column_stack([v00, v01, v11]),
         np.column_stack([v00, v10, v11]) + nv, np.column_stack([v00, v11, v01]) + nv]
    # boundary walls: directed loop of each top quad is v00->v01->v11->v10
    loops = [(v00, v01, (-1, 0)), (v01, v11, (0, 1)), (v11, v10, (1, 0)), (v10, v00, (0, -1))]
    pad = np.pad(cells, 1)
    for a, b, (di, dj) in loops:
        nb = pad[ci + 1 + di, cj + 1 + dj]
        sel = ~nb
        A, B = a[sel], b[sel]
        F.append(np.column_stack([B, A, A + nv]))
        F.append(np.column_stack([B, A + nv, B + nv]))
    return trimesh.Trimesh(V, np.vstack(F), process=True)


def to_manifold(tm):
    return m3d.Manifold(m3d.Mesh(vert_properties=np.asarray(tm.vertices, np.float32),
                                 tri_verts=np.asarray(tm.faces, np.uint32)))


def from_manifold(man):
    mm = man.to_mesh()
    # process=False: manifold output is already indexed; trimesh's tolerance-based vertex
    # merge can fuse near-coincident vertices and create non-manifold edges.
    return trimesh.Trimesh(np.asarray(mm.vert_properties)[:, :3], np.asarray(mm.tri_verts), process=False)


def _frame(n):
    n = np.asarray(n, float)
    u = np.array([1.0, 0, 0]) - n[0] * n
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    if v[1] < 0:
        v = -v
        u = -u
    return u, v


def oriented_box(center, axes, size):
    """Box with half-sizes `size` along orthonormal `axes` (3 vectors) at `center`."""
    box = m3d.Manifold.cube([2 * s for s in size], center=True)
    R = np.column_stack(axes)
    M = np.hstack([R, np.asarray(center, float)[:, None]])
    return box.transform(M.tolist())


def build_guide(bone, tumor, plan, p: GuideParams):
    n = np.asarray(plan["n"])
    d = plan["d"]
    qz = plan["q_plane"][2]
    g = p.grid_mm
    x0 = round(plan["q_plane"][0] / g) * g          # centre on the patient's sacral midline
    xs = np.arange(x0 - p.lateral_mm, x0 + p.lateral_mm + g / 2, g)
    zs = np.arange(qz - p.band_mm - 25, qz + p.band_mm + 25 + g / 2, g)
    surf, hitlists = _first_hits(bone, xs, zs)
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    S = np.stack([X, surf, Z], -1)
    h = S @ n - d

    bonehit = np.isfinite(surf)
    filled = _fill_nan(np.where(bonehit, surf, np.nan))
    ref = ndi.median_filter(filled, size=int(8 / g) | 1)
    hole = (~bonehit) | (filled < ref - 2.0)
    hole = ndi.binary_dilation(hole, iterations=int(1.0 / g))  # stay 1 mm off hole edges

    # tumor keep-out: no contact on bone within keepout of the tumor surface
    near_tumor = np.zeros_like(hole)
    tree = cKDTree(tumor.sample(60000) if tumor.area > 0 else tumor.vertices)
    if bonehit.any():
        dist, _ = tree.query(S[bonehit])
        near_tumor[bonehit] = dist < p.tumor_keepout_mm

    hh = np.where(np.isfinite(h), h, (np.stack([X, filled, Z], -1) @ n - d))
    foot = (np.abs(hh) <= p.band_mm) & np.isfinite(filled)
    foot &= ndi.binary_dilation(bonehit, iterations=2)          # only over the bone
    lab, k = ndi.label(foot)
    if k > 1:
        sizes = ndi.sum(foot, lab, range(1, k + 1))
        foot = lab == 1 + np.argmax(sizes)
    foot = ndi.binary_fill_holes(foot)
    contact = foot & bonehit & ~hole & ~near_tumor

    bottom = filled + p.clearance_mm + np.where(contact, 0.0, p.hole_relief_mm)
    bottom = np.where(contact, bottom, np.maximum(bottom, ndi.gaussian_filter(bottom, 2)))
    top = ndi.gaussian_filter(bottom, 3.0 / g) + p.base_thickness_mm
    top = np.maximum(top, bottom + 3.0)
    top_pts = np.stack([X, top, Z], -1)
    in_housing = np.abs(top_pts @ n - d) <= p.housing_half_mm
    top = top + np.where(in_housing, p.housing_height_mm, 0.0)

    cells = foot[:-1, :-1] & foot[1:, :-1] & foot[:-1, 1:] & foot[1:, 1:]
    cells = ndi.binary_opening(cells, iterations=2)
    cells = _fix_pinches(cells)
    base = heightmap_solid(xs, zs, bottom, top, cells)
    man = to_manifold(base)

    # ---- slot: exact plane, closed ends
    u, v = _frame(n)
    on_line = cells & (np.abs(hh[:-1, :-1]) < 1.0)
    cx = (xs[:-1][:, None] + 0 * zs[None, :-1])[on_line]
    cz = (zs[None, :-1] + 0 * xs[:-1][:, None])[on_line]
    cy = top[:-1, :-1][on_line]
    P = np.column_stack([cx, cy, cz])
    us = P @ u
    u0, u1 = us.min() + p.end_bridge_mm, us.max() - p.end_bridge_mm
    c_on_plane = n * d
    center = c_on_plane + u * ((u0 + u1) / 2 - c_on_plane @ u) + v * ((P @ v).mean() - c_on_plane @ v)
    half_w = (p.blade_mm + p.slot_clearance_mm) / 2
    slot = oriented_box(center, [u, v, n], [(u1 - u0) / 2, 60.0, half_w])

    kw = []
    if p.n_kwires:
        kw = choose_kwires(bone, tumor, plan, p, xs, zs, contact, hitlists, top, hh)
        for (x, z, ytop, yent) in kw:
            boss = m3d.Manifold.cylinder(ytop + p.kwire_boss_h_mm - (yent - 2.0), p.kwire_boss_r_mm,
                                         p.kwire_boss_r_mm, 48)
            boss = boss.transform([[1, 0, 0, x], [0, 0, 1, yent - 2.0], [0, -1, 0, z]])
            # clip boss so it cannot intrude below the guide underside (bone side)
            man = man + (boss - bone_shadow(bone_hm=(xs, zs, bottom), center=(x, z), r=p.kwire_boss_r_mm + 1))
        for (x, z, ytop, yent) in kw:
            hole_c = m3d.Manifold.cylinder(200.0, p.kwire_hole_mm / 2, p.kwire_hole_mm / 2, 32)
            hole_c = hole_c.transform([[1, 0, 0, x], [0, 0, 1, -100.0], [0, -1, 0, z]])
            man = man - hole_c
    man = man - slot
    mesh = from_manifold(man)
    return Guide(p, xs, zs, bottom, top, cells, contact, surf, h, (u0, u1), kw, mesh)


def bone_shadow(bone_hm, center, r):
    """Solid occupying everything BELOW the guide underside near (x,z): used to trim bosses."""
    xs, zs, bottom = bone_hm
    x, z = center
    i0, i1 = np.searchsorted(xs, [x - r, x + r])
    j0, j1 = np.searchsorted(zs, [z - r, z + r])
    sub_b = bottom[i0:i1 + 1, j0:j1 + 1]
    lo = np.full_like(sub_b, np.nanmin(sub_b) - 80.0)
    cells = np.ones((sub_b.shape[0] - 1, sub_b.shape[1] - 1), bool)
    return to_manifold(heightmap_solid(xs[i0:i1 + 1], zs[j0:j1 + 1], lo, sub_b, cells))


def choose_kwires(bone, tumor, plan, p, xs, zs, contact, hitlists, top, hh):
    """Greedy, spread-out K-wire sites that satisfy the safety rules."""
    nx, nz = len(xs), len(zs)
    tree = cKDTree(tumor.sample(60000))
    clear_of_housing = np.abs(hh) > (p.housing_half_mm + p.kwire_boss_r_mm + 1.5)
    r_cells = int(np.ceil((p.kwire_boss_r_mm + 0.5) / p.grid_mm))
    solid_contact = ndi.binary_erosion(contact, iterations=r_cells)
    cand = np.argwhere(solid_contact & clear_of_housing)
    ok = []
    for i, j in cand[:: max(1, len(cand) // 4000)]:
        hl = hitlists[i * nz + j]
        if len(hl) < 2:
            continue
        thickness = hl[0] - hl[1]   # dorsal surface to next exit (canal, foramen, or anterior cortex)
        if thickness < p.kwire_depth_mm + p.kwire_safety_mm:
            continue
        tip = np.array([xs[i], hl[0] - p.kwire_depth_mm, zs[j]])
        entry = np.array([xs[i], hl[0], zs[j]])
        seg = np.linspace(entry, tip, 8)
        dist, _ = tree.query(seg)
        if dist.min() < 10.0:
            continue
        ok.append((i, j, hh[i, j]))
    if not ok:
        return []
    ok = np.array(ok)
    pts = np.column_stack([xs[ok[:, 0].astype(int)], zs[ok[:, 1].astype(int)]])
    # preferred: preserved (cephalad, h>0) side; seed with the most lateral cephalad site
    score = np.where(ok[:, 2] > 0, 1.0, 0.6)
    chosen = [int(np.argmax(np.abs(pts[:, 0]) * score))]
    while len(chosen) < p.n_kwires:
        dmin = np.min(np.linalg.norm(pts[:, None, :] - pts[chosen][None], axis=-1), axis=1)
        chosen.append(int(np.argmax(dmin * score)))
    out = []
    for c in chosen:
        i, j = int(ok[c, 0]), int(ok[c, 1])
        out.append((float(xs[i]), float(zs[j]), float(top[i, j]), float(hitlists[i * nz + j][0])))
    return out
