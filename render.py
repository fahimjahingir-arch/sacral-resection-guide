"""Figures for the report. 3D views with PyVista (VTK, off-screen); charts with matplotlib."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

BONE, TUMOR, GUIDE, PLANE = "#e3d5b8", "#c0504d", "#3f7fbf", "#2e8b57"
pv.global_theme.font.color = "black"


def _pv(mesh):
    f = np.hstack([np.full((len(mesh.faces), 1), 3), mesh.faces]).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, float), f)


def _plotter(size=(1100, 900)):
    p = pv.Plotter(off_screen=True, window_size=size)
    p.set_background("white")
    p.enable_anti_aliasing("ssaa")
    return p


def _cam(p, center, eye_dir, dist, up=(0, 0, 1)):
    eye_dir = np.asarray(eye_dir, float)
    eye_dir /= np.linalg.norm(eye_dir)
    p.camera_position = [tuple(center + dist * eye_dir), tuple(center), up]


def _plane_mesh(plan, size=150):
    n = np.asarray(plan["n"]); d = plan["d"]
    c = n * d
    return pv.Plane(center=c, direction=n, i_size=size, j_size=size)


def scene(bone, tumor, guide_mesh, plan, path, eye=(0.45, 1.0, 0.55), title=None, show_plane=True,
          bone_opacity=1.0, lift=0.0, dist=330):
    p = _plotter()
    p.add_mesh(_pv(bone), color=BONE, opacity=bone_opacity, smooth_shading=True, specular=0.2)
    if tumor is not None:
        p.add_mesh(_pv(tumor), color=TUMOR, smooth_shading=True, opacity=0.95)
    if guide_mesh is not None:
        p.add_mesh(_pv(guide_mesh).translate((0, lift, 0)), color=GUIDE, specular=0.3)
    if show_plane:
        p.add_mesh(_plane_mesh(plan), color=PLANE, opacity=0.35)
    _cam(p, bone.bounds.mean(0), eye, dist)
    if title:
        p.add_text(title, font_size=12, position="upper_edge")
    p.screenshot(path)
    p.close()


def guide_only(gmesh, path, title):
    imgs = []
    c = gmesh.bounds.mean(0)
    for eye, lab in [((0.5, 1.0, 0.7), "top (posterior view)"),
                     ((0.4, -1.0, -0.5), "underside: bone-matched surface")]:
        p = _plotter((800, 650))
        p.add_mesh(_pv(gmesh), color=GUIDE, specular=0.3)
        _cam(p, c, eye, 200)
        p.add_text(lab, font_size=11, position="upper_edge")
        imgs.append(p.screenshot(return_img=True))
        p.close()
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=110)
    for ax, im in zip(axs, imgs):
        ax.imshow(im)
        ax.axis("off")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def section(bone, tumor, g, plan, path, dx=-12.0):
    """Sagittal cross-section through guide, bone and tumor at lateral position x."""
    x = plan["q_plane"][0] + dx            # offset from the patient midline (avoid the crest)
    fig, ax = plt.subplots(figsize=(6, 7), dpi=130)
    for mesh, col, lw, lab in [(bone, "#8a7a5c", 1.2, "bone"), (tumor, TUMOR, 1.2, "tumor"),
                               (g.mesh, GUIDE, 1.6, "guide")]:
        sec = mesh.section(plane_origin=[x, 0, 0], plane_normal=[1, 0, 0])
        if sec is None:
            continue
        first = True
        for ent in sec.discrete:
            ax.plot(ent[:, 1], ent[:, 2], color=col, lw=lw, label=lab if first else None)
            first = False
    n = np.asarray(plan["n"]); d = plan["d"]
    yc = plan["q_plane"][1]
    ys = np.array([yc - 60.0, yc + 60.0])
    ax.plot(ys, (d - n[1] * ys) / n[2], color=PLANE, ls="--", lw=1, label="planned cut")
    ax.set_aspect("equal")
    ax.set_xlabel("y  anterior → posterior (mm)")
    ax.set_ylabel("z (mm)")
    zc = plan["q_plane"][2]
    ax.set_ylim(zc - 45, zc + 45)
    ax.set_xlim(yc - 50, yc + 55)
    ax.set_title(f"Sagittal section {dx:+g} mm from midline (x = {x:.0f} mm)", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def ct_slice(OUT, plan, path):
    f_real = os.path.join(OUT, "ct_midsagittal.npz")
    f = os.path.join(OUT, "phantom_hu_midsagittal.npy")
    if os.path.exists(f_real):
        z = np.load(f_real)
        sl, ext = z["slice"], list(z["extent"])
        title = "TCIA CT + simulated chordoma, sagittal through tumor"
    elif os.path.exists(f):
        import phantom as P
        sl = np.load(f)  # (y, z)
        ext = [P.ORIGIN[1], P.ORIGIN[1] + sl.shape[0], P.ORIGIN[2], P.ORIGIN[2] + sl.shape[1]]
        title = "Synthetic CT, mid-sagittal (W 1200 / L 400)"
    else:
        return
    fig, ax = plt.subplots(figsize=(5, 6), dpi=130)
    ax.imshow(sl.T, origin="lower", cmap="gray", vmin=-200, vmax=1000, extent=ext)
    n = np.asarray(plan["n"]); d = plan["d"]
    ys = np.array([ext[0], ext[1]])
    ax.plot(ys, (d - n[1] * ys) / n[2], color=PLANE, lw=2, label="planned osteotomy")
    ax.plot(ys, (d - plan["margin_mm"] - n[1] * ys) / n[2], color=TUMOR, lw=1, ls="--",
            label="tumor tangent (margin %.0f mm)" % plan["margin_mm"])
    ax.set_xlabel("y  anterior → posterior (mm)")
    ax.set_ylabel("z  inferior → superior (mm)")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def contact_maps(guides, path):
    cmap = matplotlib.colors.ListedColormap(["white", "#f2c46d", GUIDE])
    fig, axs = plt.subplots(1, len(guides), figsize=(4.2 * len(guides), 4), dpi=130)
    for ax, (name, g) in zip(np.atleast_1d(axs), guides.items()):
        img = np.zeros(g.contact.shape)
        foot = np.zeros_like(g.contact)
        foot[:-1, :-1] |= g.cells
        img[foot] = 1
        img[g.contact & foot] = 2
        ax.imshow(img.T, origin="lower", cmap=cmap, extent=[g.xs[0], g.xs[-1], g.zs[0], g.zs[-1]],
                  vmin=0, vmax=2)
        for (x, z, *_r) in g.kwires:
            ax.plot(x, z, "o", mfc="none", mec="k", ms=8)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("z (mm)")
    fig.suptitle("Underside: blue = bone contact, amber = bridged (foramina / tumor keep-out / relief), "
                 "○ = K-wire", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def bar_charts(res, path):
    V = res["variants"]
    names = list(V)
    off = [V[k]["stability"]["offset_p95_mm"] for k in names]
    ang = [V[k]["stability"]["angle_p95_deg"] for k in names]
    area = [V[k]["contact_area_cm2"] for k in names]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), dpi=130)
    for ax, vals, lab in zip(axs, [area, off, ang], ["contact area (cm²)", "95th pct cut offset (mm)",
                                                    "95th pct cut tilt (°)"]):
        ax.bar(names, vals, color=GUIDE)
        ax.set_title(lab, fontsize=10)
        for i, v in enumerate(vals):
            ax.text(i, v, "%.2f" % v, ha="center", va="bottom", fontsize=9)
        ax.tick_params(axis="x", labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Simulated seating error (τ = %.1f mm RMS contact mismatch)" %
                 V[names[0]]["stability"]["tau_mm"], fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def all_figures(bone, tumor, guides, plan, res, OUT):
    fig = os.path.join(OUT, "figures")
    os.makedirs(fig, exist_ok=True)
    ct_slice(OUT, plan, os.path.join(fig, "01_ct_plan.png"))
    scene(bone, tumor, None, plan, os.path.join(fig, "02_plan_3d.png"), eye=(1.0, -0.12, 0.08),
          title="Lateral view: plane clears the tumor by the 10 mm margin", bone_opacity=0.6)
    gC = guides.get("C_wide_kwires") or list(guides.values())[-1]
    scene(bone, tumor, gC.mesh, plan, os.path.join(fig, "03_guide_on_bone.png"),
          title="Guide C seated on the dorsal sacrum", show_plane=False)
    scene(bone, tumor, gC.mesh, plan, os.path.join(fig, "07_exploded.png"), lift=35.0,
          title="Exploded: guide lifted 35 mm along the insertion direction")
    section(bone, tumor, gC, plan, os.path.join(fig, "08_section.png"))
    for name, g in guides.items():
        guide_only(g.mesh, os.path.join(fig, f"04_guide_{name}.png"), name)
    contact_maps(guides, os.path.join(fig, "05_contact_maps.png"))
    bar_charts(res, os.path.join(fig, "06_stability.png"))


def fixture_figure(mesh, n, u, v, path):
    """Test fixture standing on its base (planned plane horizontal), seen from the posterior side."""
    n, u, v = (np.asarray(x, float) for x in (n, u, v))
    p = _plotter((1000, 820))
    p.add_mesh(_pv(mesh), color="#d9d6cf", smooth_shading=False, specular=0.15)
    _cam(p, mesh.bounds.mean(0), 0.75 * u + 1.0 * v + 0.55 * n, 480, up=tuple(n))
    p.add_text("Test fixture: bone + anterior cradle + backplate + base", font_size=12, position="upper_edge")
    p.screenshot(path)
    p.close()
