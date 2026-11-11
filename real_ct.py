"""
Real-CT front end: DICOM series -> sacrum surface (TotalSegmentator + HU refinement)
-> simulated chordoma -> bone.stl / tumor.stl for run_all.py.

  python real_ct.py --dicom <folder> [--out outputs] [--tumor-level 0.68]

Coordinates: DICOM patient LPS (x = patient left +, y = posterior +, z = superior +),
which is the same convention the rest of the pipeline uses.

Steps
 1. DICOM -> NIfTI (SimpleITK)
 2. TotalSegmentator (--roi_subset sacrum) -> sacrum label
 3. Surface refinement: inside a 1-voxel dilation of the label, resample the CT to
    `iso` mm (cubic) and take the HU iso-surface of (HU > thr) filled + label. This puts
    the surface on the cortex edge with sub-voxel accuracy instead of on label voxels.
 4. Simulated chordoma (NO public sacral-chordoma CT exists): lobulated ellipsoid centred
    in the sacral body at a fractional level along the sacrum (0 = S1 endplate,
    1 = apex), expanding anteriorly, replacing the bone it occupies.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import SimpleITK as sitk
from scipy import ndimage as ndi

import segment as S

HERE = os.path.dirname(os.path.abspath(__file__))


def load_dicom(folder):
    r = sitk.ImageSeriesReader()
    files = r.GetGDCMSeriesFileNames(folder)
    r.SetFileNames(files)
    img = r.Execute()
    return img


def to_array_lps(img):
    """Return volume as (x, y, z) array plus spacing/origin, axis-aligned (asserted)."""
    D = np.array(img.GetDirection()).reshape(3, 3)
    if not np.allclose(np.abs(D), np.eye(3), atol=1e-3):
        raise ValueError("oblique acquisition: resample first")
    arr = sitk.GetArrayFromImage(img).transpose(2, 1, 0).astype(np.float32)  # (x,y,z)
    sp = np.array(img.GetSpacing())
    org = np.array(img.GetOrigin())
    for k in range(3):                      # flip negative axes so index increases with LPS
        if D[k, k] < 0:
            arr = np.flip(arr, k)
            org[k] = org[k] - sp[k] * (arr.shape[k] - 1)
    return arr, sp, org


def run_totalseg(nii, outdir, fast=False):
    exe = os.path.join(os.path.dirname(sys.executable), "TotalSegmentator")
    cmd = [exe, "-i", nii, "-o", outdir, "--roi_subset", "sacrum", "--device", "cpu"]
    if fast:
        cmd.append("--fast")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def refine_surface(ct, lab, sp, org, iso=0.6, thr=150.0):
    """Sub-voxel sacrum mask on an isotropic grid within the label's bounding box."""
    idx = np.argwhere(lab)
    lo = np.maximum(idx.min(0) - 3, 0)
    hi = np.minimum(idx.max(0) + 4, np.array(lab.shape))
    sl = tuple(slice(a, b) for a, b in zip(lo, hi))
    ct_c, lab_c = ct[sl], lab[sl].astype(np.float32)
    zoom = sp / iso
    ct_i = ndi.zoom(ct_c, zoom, order=3)
    lab_i = ndi.zoom(lab_c, zoom, order=1) > 0.5
    # only voxels within ~1 original voxel (in-plane) of the label may be added
    grow = ndi.binary_dilation(lab_i, iterations=max(1, int(round(sp[:2].min() / iso))))
    m = ((ct_i > thr) & grow) | lab_i
    m = ndi.binary_fill_holes(m)
    m = S.largest_component(ndi.binary_opening(m, iterations=1))
    org_i = org + lo * sp
    return m, iso, org_i, ct_i


def sacral_axis(mask, iso, org):
    """Centroid curve of the sacrum along z (superior->inferior) and its arc length."""
    zs = np.flatnonzero(mask.any(axis=(0, 1)))
    pts = []
    for k in zs:
        ij = np.argwhere(mask[:, :, k])
        if len(ij) < 20:
            continue
        pts.append([*(ij.mean(0) * iso + org[:2]), k * iso + org[2]])
    pts = np.array(pts)[::-1]                # S1 (top) first
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    return pts, s / s[-1]


def add_chordoma(mask, iso, org, level=0.68, rx=20.0, rz=18.0, r_ant=26.0, seed=0):
    """Asymmetric mass: large presacral (anterior) component, intraosseous posterior part
    that stops at ~75 % of the local AP bone thickness, so the dorsal cortex stays intact
    (the typical CT picture of sacral chordoma, and what a posterior guide must sit on)."""
    pts, f = sacral_axis(mask, iso, org)
    c = np.array([np.interp(level, f, pts[:, k]) for k in range(3)])
    # anterior = -y in LPS. Local AP thickness near the midline at this level
    k = int(round((c[2] - org[2]) / iso))
    ij = np.argwhere(mask[:, :, k])
    xs = ij[:, 0] * iso + org[0]
    mid = np.abs(xs - c[0]) < 8.0
    ys = ij[mid, 1] * iso + org[1]
    y_ant, y_post = ys.min(), ys.max()
    thick = y_post - y_ant
    c[1] = y_ant + 0.35 * thick
    r_post = 0.40 * thick
    nx, ny, nz = mask.shape
    X, Y, Z = np.meshgrid(*(np.arange(n) * iso + o for n, o in zip((nx, ny, nz), org)), indexing="ij")
    ry = np.where(Y < c[1], r_ant, r_post)
    r = np.sqrt(((X - c[0]) / rx) ** 2 + ((Y - c[1]) / ry) ** 2 + ((Z - c[2]) / rz) ** 2)
    lobes = 0.08 * np.sin(X / 6.0) * np.sin(Z / 7.0) + 0.06 * np.cos(Y / 5.0)
    tumor = r < (1.0 + lobes)
    return tumor, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dicom", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "outputs"))
    ap.add_argument("--work", default=None, help="folder for NIfTI + TotalSegmentator output")
    ap.add_argument("--tumor-level", type=float, default=0.68)
    ap.add_argument("--iso", type=float, default=0.6)
    a = ap.parse_args()
    work = a.work or os.path.join(a.out, "ct_work")
    os.makedirs(work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    img = load_dicom(a.dicom)
    nii = os.path.join(work, "ct.nii.gz")
    sitk.WriteImage(img, nii)
    ct, sp, org = to_array_lps(img)
    print("CT", ct.shape, "spacing", sp, "origin", org)
    segdir = os.path.join(work, "seg")
    if not os.path.exists(os.path.join(segdir, "sacrum.nii.gz")):
        run_totalseg(nii, segdir)
    lab_img = sitk.ReadImage(os.path.join(segdir, "sacrum.nii.gz"))
    lab, _, _ = to_array_lps(lab_img)
    lab = lab > 0.5
    print("sacrum label voxels", int(lab.sum()), "volume %.1f mL" % (lab.sum() * sp.prod() / 1000))

    m, iso, org_i, ct_i = refine_surface(ct, lab, sp, org, a.iso)
    tumor, tc = add_chordoma(m, iso, org_i, a.tumor_level)
    bone_m = m & ~tumor
    bone = S.mask_to_mesh(bone_m, iso, org_i, smooth_sigma=1.0)
    # keep only sizeable bone pieces (tumor can split the sacrum)
    parts = [p for p in bone.split(only_watertight=True) if p.volume > 500]
    import trimesh
    bone = trimesh.util.concatenate(parts)
    tum = S.mask_to_mesh(tumor, iso, org_i, smooth_sigma=1.0)
    bone.export(os.path.join(a.out, "bone.stl"))
    tum.export(os.path.join(a.out, "tumor.stl"))

    # mid-sagittal slice for the report figure (through the tumor centre)
    ix = int(round((tc[0] - org_i[0]) / iso))
    sl = ct_i[ix].copy()
    sl[tumor[ix]] = np.clip(sl[tumor[ix]], None, 40) * 0 + 35   # paint simulated tumor ~35 HU
    np.savez(os.path.join(a.out, "ct_midsagittal.npz"), slice=sl,
             extent=[org_i[1], org_i[1] + sl.shape[0] * iso, org_i[2], org_i[2] + sl.shape[1] * iso])
    info = {"source": "TCIA Pelvic-Reference-Data, Pelvic-Ref-003 (CC BY 3.0, doi:10.7937/TCIA.2019.woskq5oo)",
            "ct_shape": list(ct.shape), "ct_spacing_mm": sp.tolist(),
            "sacrum_label_volume_ml": float(lab.sum() * sp.prod() / 1000),
            "refined_iso_mm": iso, "bone_pieces": len(parts), "bone_watertight": bool(bone.is_watertight),
            "tumor": {"simulated": True, "level_fraction": a.tumor_level, "center_mm": tc.tolist(),
                      "volume_ml": float(tumor.sum() * iso ** 3 / 1000),
                      "bone_destroyed_ml": float((m & tumor).sum() * iso ** 3 / 1000)}}
    json.dump(info, open(os.path.join(a.out, "ct_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
