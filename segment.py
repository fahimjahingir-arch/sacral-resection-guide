"""
Threshold + morphology segmentation (the scripted equivalent of Slicer's
Threshold -> Islands -> Smoothing workflow) and surface extraction.
"""
import numpy as np
import trimesh
from scipy import ndimage as ndi
from skimage import measure


def largest_component(mask, seed=None):
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    if seed is not None and lab[seed] > 0:
        return lab == lab[seed]
    sizes = ndi.sum(mask, lab, range(1, n + 1))
    return lab == (1 + int(np.argmax(sizes)))


def segment_bone(hu, threshold=180.0, min_voxels=2000):
    m = hu > threshold
    m = ndi.binary_closing(m, iterations=1)
    m = ndi.binary_fill_holes(m)          # fills cancellous gaps, keeps open canal/foramina
    m = ndi.binary_opening(m, iterations=1)
    # keep every sizeable island (a destructive tumor can split the sacrum in two)
    lab, n = ndi.label(m)
    sizes = ndi.sum(m, lab, range(1, n + 1))
    keep = 1 + np.flatnonzero(sizes >= min_voxels)
    return np.isin(lab, keep)


def segment_tumor(hu, seed_ijk, lo=5.0, hi=45.0):
    """Seeded region growing in the soft-tissue window (chordoma ~20-40 HU on CT)."""
    sm = ndi.median_filter(hu, size=3)
    m = (sm > lo) & (sm < hi)
    m = ndi.binary_opening(m, iterations=2)   # break thin bridges to adjacent muscle
    m = largest_component(m, seed=tuple(seed_ijk))
    m = ndi.binary_closing(m, iterations=2)
    return ndi.binary_fill_holes(m)


def mask_to_mesh(mask, spacing, origin, smooth_sigma=0.8):
    vol = ndi.gaussian_filter(mask.astype(np.float32), smooth_sigma)
    verts, faces, _, _ = measure.marching_cubes(np.pad(vol, 1), level=0.5, spacing=(spacing,) * 3)
    verts = verts - spacing + origin
    mesh = trimesh.Trimesh(verts, faces, process=True)
    mesh.fix_normals()
    if mesh.volume < 0:
        mesh.invert()
    return mesh


def dice(a, b):
    return 2.0 * np.logical_and(a, b).sum() / (a.sum() + b.sum())
