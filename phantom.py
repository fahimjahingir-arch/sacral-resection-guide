"""
Synthetic sacral CT phantom with a simulated chordoma.

PURPOSE: a known-truth test volume so the pipeline can be run and checked end to end
without patient data. It is NOT anatomy-accurate. For the real project, segment a
TCIA pelvic CT in 3D Slicer, export bone.stl + tumor.stl, and run
`run_all.py --bone bone.stl --tumor tumor.stl` instead.

Coordinate frame (mm): x = patient left(-)/right(+), y = anterior(-)/posterior(+),
z = inferior(-)/superior(+). The sacral base (S1 endplate) is at z ~ 110, apex at z ~ 0.
"""
import numpy as np
from scipy import ndimage as ndi

SPACING = 1.0  # mm, isotropic
ORIGIN = np.array([-75.0, -75.0, -10.0])
SHAPE = (150, 130, 140)  # x, y, z voxels

HU_FAT, HU_MUSCLE, HU_CANCELLOUS, HU_CORTEX, HU_TUMOR = -100, 60, 250, 1200, 30

# Foramen levels (z, mm) for S1..S4 -> used later to name the cut level
FORAMEN_Z = {"S1": 95.0, "S2": 75.0, "S3": 55.0, "S4": 38.0}


def centerline_y(z):
    """Anterior concavity: mid-sacrum sits more posterior than base and apex."""
    return 8.0 - 0.006 * (z - 60.0) ** 2


def half_width(z):
    return 8.0 + 47.0 * np.clip(z / 110.0, 0, 1) ** 1.3


def ap_thickness(z):
    return 8.0 + 37.0 * np.clip(z / 110.0, 0, 1) ** 1.5


def grid():
    ix, iy, iz = [np.arange(n, dtype=np.float32) for n in SHAPE]
    X, Y, Z = np.meshgrid(ix * SPACING + ORIGIN[0], iy * SPACING + ORIGIN[1],
                          iz * SPACING + ORIGIN[2], indexing="ij")
    return X, Y, Z


def make_phantom(seed=0, tumor_center=(4.0, None, 31.0), tumor_radii=(24.0, 28.0, 19.0)):
    rng = np.random.default_rng(seed)
    X, Y, Z = grid()
    yc = centerline_y(Z)
    w = half_width(Z)
    t = ap_thickness(Z)
    # AP thickness reduces toward the lateral masses / ala
    t_loc = t * (1.0 - 0.45 * np.clip(np.abs(X) / w, 0, 1) ** 2)
    in_z = (Z > 0) & (Z < 110)
    body = in_z & (np.abs(X) < w) & (np.abs(Y - yc) < t_loc / 2)

    dorsal = yc + t_loc / 2
    # median sacral crest with spinous tubercles
    tub = 3.0 + 3.0 * np.cos((Z - 95.0) * 2 * np.pi / 20.0) ** 2
    crest = (Z > 30) & (Z < 104) & (np.abs(X) < 3.5) & (Y >= dorsal - 1) & (Y < dorsal + tub)
    # intermediate crests (articular tubercles), smaller
    icrest = (Z > 30) & (Z < 100) & (np.abs(np.abs(X) - 14) < 2.5) & (Y >= dorsal - 1) & (Y < dorsal + 2.5)
    bone = body | crest | icrest

    # sacral canal (open at the base, closes toward the hiatus)
    r_can = np.clip(3.0 + 6.0 * (Z - 20) / 90.0, 0, None)
    canal = (Z > 20) & (Z < 115) & (((X / 1.6) ** 2 + (Y - (yc + 0.18 * t)) ** 2) < r_can ** 2)
    bone &= ~canal

    # 4 pairs of foramina, through the full AP thickness
    foramina = np.zeros_like(bone)
    for i, fz in enumerate(FORAMEN_Z.values()):
        r = 5.0 - 0.6 * i
        fx = 0.42 * half_width(np.array(fz))
        for sx in (-1, 1):
            foramina |= ((X - sx * fx) ** 2 + (Z - fz) ** 2) < r ** 2
    bone &= ~foramina

    # ---- simulated chordoma: lobulated, anterior, midline, destroys bone it occupies
    cx, cy, cz = tumor_center
    if cy is None:
        cy = float(centerline_y(np.array(cz))) - 6.0
    rx, ry, rz = tumor_radii
    r = np.sqrt(((X - cx) / rx) ** 2 + ((Y - cy) / ry) ** 2 + ((Z - cz) / rz) ** 2)
    lobes = 0.08 * np.sin(X / 6.0) * np.sin(Z / 7.0) + 0.06 * np.cos(Y / 5.0)
    tumor = r < (1.0 + lobes)

    # ---- assemble HU
    hu = np.full(SHAPE, HU_FAT, np.float32)
    muscle = (Y > dorsal) & (Y < dorsal + 22) & (np.abs(X) < w + 10) & (Z > -5) & (Z < 125)
    hu[muscle] = HU_MUSCLE
    cortex = bone & ~ndi.binary_erosion(bone, iterations=2)
    hu[bone] = HU_CANCELLOUS
    hu[cortex] = HU_CORTEX
    hu[tumor] = HU_TUMOR
    hu = ndi.gaussian_filter(hu, 0.7)  # scanner PSF
    hu += rng.normal(0, 15, SHAPE).astype(np.float32)  # quantum noise

    truth = {"bone": bone & ~tumor, "tumor": tumor}
    return hu, truth


def voxel_to_world(ijk):
    return np.asarray(ijk) * SPACING + ORIGIN
