# Bench-test protocols

These three tests turn the design into a validated design. Each has a printable part, a data
sheet, and an analysis script in this folder. Run the pipeline first (`python ../run_all.py`) so
`outputs/` holds the bone model, guides and `results.json`.

Run the cheap pilot first (1 operator, 3 placements per variant) to check the method before you
collect the full data set.

---

## TEST 1: Cut accuracy (does the cut land where it was planned?)

**Prints**
| Part | File | Material / process | Notes |
|---|---|---|---|
| Sacrum test fixture | `outputs/test_fixture.stl` (~500 cm³ solid, ~150–200 g at 15 % infill) | PLA or PETG FDM, 0.2 mm layers, 15 % gyroid infill, 4 walls | Print **backplate-down**. The base face (the other flat face) is the measuring datum: check it with a straightedge, and sand it flat if it warped |
| Guide C | `outputs/guide_C_wide_kwires.stl` | Same process you'd propose clinically (PA12 SLS through a print service, or SLA resin) | Orient the slot vertically so the slot walls aren't stair-stepped |

**Equipment:** oscillating saw with a 1.27 mm blade (a hobby multi-tool works for printed bone), 1.6 mm
K-wires or 1.6 mm drill-rod pins, a 1.27 mm steel shim, calipers or a height gauge, a flat
table (glass or granite), and a phone with a free 3D-scanning app (optional).

**Procedure (n ≥ 5 specimens; print 5 fixtures)**
1. Lay the fixture on its **backplate**. The dorsal sacrum now faces up, as in a prone patient. Seat the guide by pressing
   straight down (the insertion direction), then insert the 3 K-wires.
2. Cut through the slot. Finish the lateral few mm freehand in the same plane, then remove the guide.
3. **Measure the cut.**
   - *With a scanner:* scan the preserved (cranial) piece and run `python analyze_cut.py --scan piece.stl`.
     If the scan's coordinate frame is far from the model, align them roughly first (3 matching points in MeshLab).
   - *Without a scanner:* stand the fixture on its **base**. The cut face is parallel to the table when the plan is perfect. Measure the height of the
     cut face above the table at 4 points (left and right, front and back) with the height gauge.
     Offset = mean height − `table_to_planned_plane_mm` (in `outputs/fixture_info.json`). Tilt comes from the
     height differences divided by the point spacing.
4. Record offset (mm), angle (°) and cut-face flatness.

**What to report:** mean ± SD offset, mean ± SD angle, worst case, and margin remaining = 10 mm + offset
(negative offset means the cut moved toward the tumor).

The analysis script was checked on synthetic cuts with known error (`python analyze_cut.py --selftest`,
results in `outputs/selftest_analyze_cut.json`).

---

## TEST 2: Placement repeatability (A vs B vs C, multiple users)

**Prints:** 1 test fixture (re-used, no cutting) + guides A, B, C.

**Procedure**
1. Recruit 3–5 people. A mix of engineers and non-engineers is fine. Brief each person for 1 minute: "press the guide
   down until it stops rocking". Use no further coaching.
2. For each variant, in randomized order: each person seats the guide 10 times. Lift the guide fully off between trials.
3. For each placement, put the 1.27 mm shim in the slot (it rests on the caudal wall). Measure shim height at the
   4 marked points with the height gauge. Variant C: insert the K-wires before measuring, and pull them between trials.
4. Fill `repeatability_datasheet.csv` (`python analyze_repeatability.py --make-template`):
   offset_mm = mean height − `expected_blade_lower_face_mm`; angle_deg = plane tilt.
5. `python analyze_repeatability.py repeatability_datasheet.csv`

**What to report:** within-operator SD (repeatability), between-operator SD (reproducibility),
95 % limits of agreement per variant, and Levene's test for a difference in variance between variants.
Compare with the simulated prediction in `outputs/results.json → variants.*.stability`.

---

## TEST 3: Sterilization dimensional stability

**Prints:** `outputs/autoclave_coupon.stl`, 5 coupons per material
(e.g. PA12 SLS and a biocompatible SLA resin; print-service quotes are usually cheap for small parts).

**Procedure**
1. Measure every feature in `autoclave_datasheet.csv` (`python autoclave_coupon.py --make` creates it)
   with calipers / pin gauges. Use 1.40 mm and 1.80 mm gauge pins if you can get them.
   For rib bow, put the rib on the table and measure the gap at the tip with feeler gauges.
2. Sterilize with steam: a standard gravity or pre-vacuum cycle (e.g. 121 °C/30 min or 134 °C/3–4 min per the
   material datasheet). A university biology or microbiology lab autoclave is the usual way to get access.
   Ask the lab manager. **Use only a cycle the material datasheet allows.**
3. Re-measure after cooling to room temperature for 24 h.
4. `python autoclave_coupon.py autoclave_datasheet.csv`

**What to report:** % change per feature and whether the 1.27 mm blade still fits the slot after sterilization
(the script flags clearance < 0.05 mm). Report warpage of the thin rib separately.

---

## Pre-registered pass criteria (write these down BEFORE testing)
| Metric | Target | Rationale |
|---|---|---|
| Cut offset, worst case | ≤ 2 mm toward tumor | Leaves ≥ 8 mm of the 10 mm planned margin |
| Cut angle, mean | ≤ 3° | Engineering target; at 40 mm from the guide, 3° ≈ 2 mm of margin |
| Placement SD (offset) | ≤ 1 mm | Engineering target; keeps 95 % of placements within 2 mm |
| Slot width after autoclave | ≥ 1.32 mm | Blade (1.27 mm) must not bind |
