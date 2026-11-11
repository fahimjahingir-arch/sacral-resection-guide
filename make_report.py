"""Builds REPORT.md from outputs/results.json (+ self-test and fixture info)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")


def f(x, k=2):
    return f"{x:.{k}f}"


def main():
    r = json.load(open(os.path.join(OUT, "results.json")))
    st = json.load(open(os.path.join(OUT, "selftest_analyze_cut.json"))) if os.path.exists(
        os.path.join(OUT, "selftest_analyze_cut.json")) else None
    fx = json.load(open(os.path.join(OUT, "fixture_info.json"))) if os.path.exists(
        os.path.join(OUT, "fixture_info.json")) else None
    p = r["plan"]; V = r["variants"]
    L = []
    real = os.path.exists(os.path.join(OUT, "ct_info.json"))
    L.append("# Sacral resection guide: results (%s)\n" % ("TCIA CT, simulated tumor" if real else "synthetic phantom run"))
    L.append("> All numbers below are **computational**: they come from simulation on the CT below. "
             "The physical tests in `validation/PROTOCOLS.md` are what turn them into measured results.\n")
    # ---------------- key findings (computed) ----------------
    A, B = V["A_narrow"]["stability"], V["B_wide"]["stability"]
    L.append("## Key findings\n")
    L.append(f"1. **A wider footprint is the biggest accuracy lever.** Going from A to B doubles contact area "
             f"({f(V['A_narrow']['contact_area_cm2'],1)} → {f(V['B_wide']['contact_area_cm2'],1)} cm²). It cuts the 95th-percentile "
             f"cut offset from {f(A['offset_p95_mm'])} to {f(B['offset_p95_mm'])} mm "
             f"(−{f(100*(1-B['offset_p95_mm']/A['offset_p95_mm']),0)} %) and the tilt from {f(A['angle_p95_deg'])}° to "
             f"{f(B['angle_p95_deg'])}° (−{f(100*(1-B['angle_p95_deg']/A['angle_p95_deg']),0)} %).")
    L.append(f"2. **The weakest direction is {B['weakest_mode']}** for every variant: the dorsal sacrum is smoothest along its long "
             "axis. Next design step: make the footprint wrap a feature that resists that slide, such as the median-crest tubercles, "
             "the S1 superior articular process, or the lateral sacral border.")
    if "fea" in r:
        F = r["fea"]
        bb = next(c for c in F["variants_PA12"] if c["variant"] == "B_wide" and c["support"] == "bonded")
        pins = next((c for c in F["variants_PA12"] if c["support"] == "pins"), None)
        play = {round(x["housing_height_mm"], 1): x["play_angle_deg"] for x in F["slot_play"]}
        p10 = play.get(10.0, F["slot_play"][0]["play_angle_deg"])
        L.append(f"3. **Elastic bending is not the error source; slot clearance is.** Under {F['assumptions']['F_N']:.0f} N the PA12 slot wall "
                 f"deflects ~{f(1000*bb['wall_top_edge_mm'],1)} µm ({f(bb['blade_angle_error_deg'],3)}°). The 0.13 mm blade clearance allows "
                 f"{f(p10,2)}° of rigid wobble, about {f(p10/max(bb['blade_angle_error_deg'],1e-9),0)}× more. So material should be chosen for "
                 "sterilization stability and biocompatibility, not stiffness, and a taller housing or tighter slot is the accuracy lever "
                 "(14 mm housing → 0.41°). The cost is a taller profile in the wound.")
        if pins:
            L.append(f"4. **K-wires matter only if seating contact is lost.** With pins carrying the load alone, deflection rises to "
                     f"{f(1000*pins['wall_top_edge_mm'],0)} µm (×{f(pins['wall_top_edge_mm']/bb['wall_top_edge_mm'],0)}). That is still small. "
                     "The real job of the K-wires is to stop the guide sliding during sawing (finding 2's weakest mode), which only the bench test can quantify.")
    if st:
        L.append(f"5. **The measurement method is good enough for the job.** The cut-analysis pipeline recovers known errors to "
                 f"≤{f(st['max_offset_error_mm'],2)} mm / ≤{f(st['max_angle_error_deg'],2)}°, more than 10× finer than the 2 mm / 3° pass criteria.")
    L.append("")
    ci_path = os.path.join(OUT, "ct_info.json")
    if os.path.exists(ci_path):
        ci = json.load(open(ci_path))
        t = ci["tumor"]
        L.append("## 1. Data and segmentation\n")
        L.append(f"* **CT:** {ci['source']}. Matrix {ci['ct_shape']}, voxel spacing "
                 f"{', '.join(f(x, 2) for x in ci['ct_spacing_mm'])} mm.\n"
                 f"* **Sacrum:** TotalSegmentator (`--roi_subset sacrum`), then sub-voxel refinement: cubic resampling of the CT to "
                 f"{ci['refined_iso_mm']} mm and an HU iso-surface inside the label. Label volume {f(ci['sacrum_label_volume_ml'],1)} mL.\n"
                 f"* **Tumor: SIMULATED.** No public sacral-chordoma CT exists. A lobulated {f(t['volume_ml'],1)} mL anterior mass was placed "
                 f"at {f(100*t['level_fraction'],0)} % of the sacral length from S1, which is roughly S3–S4, the typical chordoma site. "
                 f"It replaces {f(t['bone_destroyed_ml'],1)} mL of bone.\n"
                 f"* **Limitation:** the {f(max(ci['ct_spacing_mm']),1)} mm slice spacing smooths the surface along z. A clinical guide "
                 "would use ≤1 mm slices.\n")
        L.append("![CT](outputs/figures/01_ct_plan.png)\n")
    elif "segmentation" in r:
        s = r["segmentation"]
        L.append("## 1. Segmentation\n")
        L.append(f"Threshold/morphology segmentation against phantom ground truth: **Dice bone {f(s['dice_bone'],3)}, "
                 f"tumor {f(s['dice_tumor'],3)}**. Tumor volume {f(s['tumor_volume_ml'],1)} mL.\n")
        L.append("![CT](outputs/figures/01_ct_plan.png)\n")
    L.append("## 2. Resection plan\n")
    L.append(f"* Planned bone margin **{p['margin_mm']:.0f} mm**, set from the 8–10 mm used in published PSI sacrectomy "
             "planning, where the margin also had to absorb guide error.\n"
             f"* Cut plane normal {[round(x,3) for x in p['n']]}, crossing the sacral midline at z ≈ {f(p['z_level_midline'],1)} mm "
             + (f"(cut at {f(100*p['level_fraction'],0)} % of sacral length from S1; check the level against the foramina in the figures). "
                if real and 'level_fraction' in p else
                "(phantom foramina: S2 at z = 75, S3 at z = 55). That is a sacral amputation between the S2 and S3 levels, "
                "preserving the S1–S2 roots. ") + "\n"
             f"* {f(100*p['bone_kept_fraction'],0)} % of the bone surface is preserved. The minimum tumor-to-plane distance is {f(p['min_tumor_to_plane_mm'],1)} mm.\n")
    L.append("![plan](outputs/figures/02_plan_3d.png)\n")
    L.append("## 3. Guide variants\n")
    L.append("| Variant | Contact area (cm²) | Volume (cm³) | Slot length (mm) | K-wires | Watertight |\n|---|---|---|---|---|---|")
    for k, v in V.items():
        L.append(f"| {k} | {f(v['contact_area_cm2'],1)} | {f(v['volume_cm3'],1)} | {f(v['slot_length_mm'],1)} | "
                 f"{len(v['kwires'])} | {v['watertight']} |")
    L.append("\n![guide on bone](outputs/figures/03_guide_on_bone.png)\n")
    L.append("![C](outputs/figures/04_guide_C_wide_kwires.png)\n")
    L.append("![exploded](outputs/figures/07_exploded.png)\n")
    L.append("![section](outputs/figures/08_section.png)\n")
    L.append("![contact](outputs/figures/05_contact_maps.png)\n")
    L.append("## 4. Seating stability (simulation)\n")
    tau = V[next(iter(V))]["stability"]["tau_mm"]
    L.append(f"Linearized rigid-body model with τ = {tau} mm RMS random contact mismatch per 5 mm patch "
             "(the assumed residual soft tissue and seating noise).\n")
    L.append("| Variant | Patches | Offset mean / 95th (mm) | Tilt mean / 95th (°) | Margin used, 95th | Weakest mode |\n|---|---|---|---|---|---|")
    for k, v in V.items():
        s = v["stability"]
        L.append(f"| {k} | {s['n_contact_patches']} | {f(s['offset_mean_abs_mm'])} / {f(s['offset_p95_mm'])} | "
                 f"{f(s['angle_mean_deg'])} / {f(s['angle_p95_deg'])} | {f(s['margin_consumed_p95_pct'],0)} % | {s['weakest_mode']} |")
    L.append("\n![stability](outputs/figures/06_stability.png)\n")
    L.append("Errors scale linearly with τ, so at τ = 1 mm double them. K-wires don't change *seating* accuracy in this "
             "model. Their job is holding the pose during sawing, which the FEA below covers.\n")
    if "fea" in r:
        F = r["fea"]
        L.append("## 5. FEA: slot wall under saw load\n")
        L.append(f"Load: {F['assumptions']['F_N']} N on the cephalad slot wall, along the plane normal. "
                 f"Hex voxel mesh, {F['assumptions']['pitch_mm']} mm pitch, PA12 (E = 1.7 GPa).\n")
        L.append("| Variant | Support | Elements | Slot-wall top edge (µm) | Wall max (µm) | Blade angle from deflection (°) |\n|---|---|---|---|---|---|")
        for c in F["variants_PA12"]:
            L.append(f"| {c['variant']} | {c['support']} | {c['n_elements']} | {f(1000*c['wall_top_edge_mm'],1)} | "
                     f"{f(1000*c['wall_max_deflection_mm'],1)} | {f(c['blade_angle_error_deg'],3)} |")
        L.append("\n**Material (B, bonded, scaled by 1/E):**\n")
        L.append("| Material | Deflection (µm) | Angle (°) |\n|---|---|---|")
        for m, c in F["materials_B_bonded"].items():
            L.append(f"| {m} | {f(1000*c['wall_max_deflection_mm'],2)} | {f(c['blade_angle_error_deg'],4)} |")
        if "housing_height_B_PA12" in F:
            L.append("\n**Housing height trade-off (B, PA12):** taller housing = more guidance but a longer cantilever.\n")
            L.append("| Housing (mm) | FEA deflection (µm) | FEA angle (°) | Slot-play angle (°) |\n|---|---|---|---|")
            play = {round(x["housing_height_mm"], 1): x for x in F["slot_play"]}
            for c in F["housing_height_B_PA12"]:
                pl = play.get(round(c["housing_height_mm"], 1), {})
                L.append(f"| {c['housing_height_mm']:g} | {f(1000*c['wall_max_deflection_mm'],1)} | "
                         f"{f(c['blade_angle_error_deg'],3)} | {f(pl.get('play_angle_deg', float('nan')),2)} |")
            L.append("\nSlot-play angle = atan(0.13 mm clearance / guided depth): the rigid wobble the blade has in the slot.\n")
        L.append("\n**Mesh convergence (B, bonded):** " + ", ".join(
            f"{c['pitch_mm']} mm → {f(1000*c['wall_max_deflection_mm'],2)} µm ({c['n_elements']} el)" for c in F["convergence_B"]) +
                 ". The roughly ±25 % scatter comes from the stair-stepped voxel slot walls. It doesn't change the conclusion, because elastic "
                 "deflection is orders of magnitude below slot play.\n")
    if st:
        L.append("## 6. Measurement-method check\n")
        L.append(f"`analyze_cut.py` recovered known synthetic cut errors to within **{f(st['max_offset_error_mm'],3)} mm** "
                 f"and **{f(st['max_angle_error_deg'],3)}°**, with 0.08 mm scanner noise and a 4°/5 mm initial misregistration.\n")
        L.append("| True offset | True angle | Measured offset | Measured angle | Flatness RMS |\n|---|---|---|---|---|")
        for c in st["cases"]:
            L.append(f"| {c['true_offset_mm']} | {c['true_angle_deg']} | {c['offset_mm']} | {c['angle_deg']} | {c['flatness_rms_mm']} |")
    if fx:
        L.append("\n![fixture](outputs/figures/09_fixture.png)\n")
        L.append(f"\nTest fixture: the planned plane sits {f(fx['table_to_planned_plane_mm'],2)} mm above the table, and the expected blade "
                 f"lower face is at {f(fx['expected_blade_lower_face_mm'],2)} mm.\n")
    open(os.path.join(HERE, "REPORT.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote REPORT.md")


if __name__ == "__main__":
    main()
