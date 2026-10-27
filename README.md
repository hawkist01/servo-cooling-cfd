# Servo Motor Conjugate Heat Transfer — Air-Cooled Finned Housing

Steady-state conjugate heat transfer (CHT) of a robot-joint servo motor: 80 W of
copper loss generated inside the winding, conducted out through a closed
aluminium housing with 12 axial fins, and carried away by 3 m/s of air flowing
along the motor axis. The point of the study is not the final temperature — it is
the decomposition of *where* the temperature rise happens.

| | |
|---|---|
| Solver | SimScale, Conjugate Heat Transfer (fluid + 2 solid regions) |
| Turbulence | k-ω SST, steady RANS |
| Coolant | Air @ 25 °C, 3 m/s axial, 0.323 kg/s |
| Winding | copper Ø52 × 112 mm, 80 W volumetric source (336 473 W/m³) |
| Housing | aluminium Ø64 × 116 mm, closed ends, 12 axial fins 10 × 1.5 mm |
| Air domain | 600 × 300 × 300 mm, inlet 150 mm upstream |
| Mesh | 683 624 cells / 204 600 nodes, 3 boundary layers |
| Runs | 2 (see §8 of the report for why) |

## Headline result

With the winding peak at **79.9 °C** and the wetted surface **54.2 K** above the
inlet air, **54.0 K of the rise is the air film**, **0.37 K is the metal**, and
**0.25 K is the air's own bulk rise**. The case is air-side limited: no amount of
copper or aluminium moves the answer.

The pre-run lumped model predicted an 87.9 °C hot spot, so the CFD came in
**12.6 % low on superheat**. The effective coefficient is **h = 26.0 W/m²K**
against the flat-plate anchor's 22.5 (+15.9 %), and the thermal resistance
**0.678 K/W** against 0.784 (−13.6 %). The correlation is conservative here, and
the reason is physical: 12 short developing channels fed through a 180° turn
transfer better than a flat isothermal plate does.

The solid is very nearly isothermal. The whole housing spans **2.8 K** from its
coldest fin tip to its hottest patch, and the winding spans **0.87 K**. With
Bi ≈ 1.1 × 10⁻³ the aluminium is lumped by any standard, and the fin efficiency
is 0.995, so the fins are as hot as the barrel they grow from.

## Why the energy balance looks wrong, and is not

The area-average outlet temperature implies 101.1 W leaving an 80 W source, a
+26.3 % "imbalance". It is not an imbalance. An area average weights the slow warm
wake behind the motor as heavily as the fast cold core; the mass-weighted rise is
0.2466 K by definition. The gap is the velocity–temperature covariance at the
outlet. The checks that *are* exact: inlet T to six figures, mass conserved to
0.0012 % between inlet and outlet, outlet lateral velocity ~3 × 10⁻⁵ m/s. See
`docs/report.md` §4.1.

## Figures

The part, to scale — a closed housing with 12 axial fins around the winding:

![Geometry to scale](geometry/geometry_preview.png)

Where the temperature actually goes, and how the pre-run prediction compares:

| Thermal path | Prediction vs CFD |
|---|---|
| ![Thermal path](figures/thermal_path.png) | ![Prediction vs CFD](figures/validation_prediction.png) |

The energy-balance artefact and the convergence history:

| Energy balance | Convergence |
|---|---|
| ![Energy balance](figures/energy_balance.png) | ![Convergence](figures/convergence.png) |

Full write-up, including the uncertainty budget and what each omission would cost,
is in [`docs/report.md`](docs/report.md). A single-file summary is
[`docs/dashboard.html`](docs/dashboard.html).

## Layout

```
geometry/
  make_geometry.py         three-solid STEP + STLs + design_point.json
  preview.py               two-panel scale drawing (geometry_preview.png)
  design_point.json        the prediction, fixed before the mesh existed
analysis/
  derive.py                measurements.json + design_point.json -> summary.json
  make_figures.py          the four figures
  fill_report.py           template + summary.json -> docs/report.md
  make_dashboard.py        single-file HTML summary
results/
  measurements.json        converged monitors, solid field statistics, provenance
  convergence.json         downsampled monitors and residuals
  summary.json             generated derived quantities
docs/
  report_template.md       report prose with numeric placeholders (source of truth)
  report.md                generated report — do not hand-edit
  dashboard.html           generated summary page
build.py                   derive -> figures -> report -> dashboard
figures/                   generated PNGs
LICENSE                    MIT
```

## Reproduce

```
python geometry/make_geometry.py     # STEP + STLs + design_point.json
python geometry/preview.py           # geometry/geometry_preview.png
python build.py                      # derive -> figures -> report -> dashboard
```

`build.py` regenerates `results/summary.json`, `figures/*.png`, the report and the
dashboard from `results/measurements.json`, and stops on the first failure — so a
stale report cannot sit next to fresh numbers. Every number in the report is a
`{{placeholder}}` resolved from `summary.json`; a missing key is a hard error, not
a blank.

`derive.py` does not re-derive CFD. The fluid values are area-average result
controls read from the SimScale results API at the final iteration. The solid
values come from the solution archive itself: `result.zip` holds OpenFOAM ASCII
fields per region, and `1000/region2/T.gz` / `1000/region3/T.gz` were decompressed
and parsed. `region2` has 66 boundary patches and `region3` has 3 — exactly the
housing's and the winding's CAD face counts, which is what proves the mapping.
Those are cell-centre values, so the solid averages are volume means.

## Verification practice

Seven problems were caught before they reached a published number, and none of them
by "looks reasonable". The three most transferable:

1. **The STEP was 1000× too small.** The OCC writer declared `SI_UNIT(.MILLI.)`
   while writing metre-valued coordinates. gmsh round-trips agree with themselves,
   so only parsing the *declared unit* against the *coordinate range* catches it.
   `make_geometry.py` now asserts this on every build.

2. **A selection test silently destroyed the outlet monitor.** Clicking a solid's
   parent row left the outlet temperature control holding 66 housing faces, with no
   UI complaint and a plausible name. It was caught by reading the persisted spec
   and counting entities before the run. Selection experiments belong on a scratch
   control; entity counts belong in the pre-run check.

3. **The solid result controls produced nothing.** Two `AREA_INTEGRAL` controls
   saved, validated and were accepted by the run, and returned no plots at all —
   the manifest held exactly 18 PLOTs, 9 fields × 2 fluid controls. A result
   control that saves is not the same as one that produces data. The solid fields
   were read from the archive instead of spending a third run on a guess.

The general lesson running through all seven: **verify the persisted spec, not the
screen.** See `docs/report.md` §6.

## Honest limitations

- **Perfect contact at 66 interfaces.** Real winding-to-housing heat transfer passes
  through slot liner, varnish and an interference-fit joint. A 1 mm air gap over the
  0.0225 m² interface adds ~1.7 K/W — more than twice the entire convective
  resistance. The winding temperature here is a **lower bound**, possibly by tens of
  kelvin, and this is the single largest idealisation in the model.
- **One mesh, no independence study.** The cold-plate study in this portfolio has
  three levels; this one has a single 683 624-cell mesh for runtime reasons. A
  10–20 % mesh-induced uncertainty on h is plausible and is not quantified. The
  agreement in §4.3 is consistent-within-uncertainty, not a measurement of the
  correlation's error.
- **Buoyancy not modelled** (no buoyancy switch in this solver version). Ri = 0.024
  at the measured superheat, so forced convection dominates, but the omission is
  real between the fins at low speed.
- **Radiation not modelled.** At the measured surface temperature the black-body
  flux is 3 % of the convective flux for bare aluminium (ε ≈ 0.1) but 24 % for an
  anodised or painted housing. The model implicitly assumes unfinished metal.
- **The correlation is used twice.** The prediction and the CFD are compared through
  a coefficient the prediction assumed, so the comparison tests the geometry and the
  resistance network, not the correlation.

## SimScale project

https://www.simscale.com/workbench/?pid=8680246712396085531

## Licence

MIT — see [LICENSE](LICENSE).
