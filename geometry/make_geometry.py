"""
Robot joint servo motor -- conjugate heat transfer geometry.

Three closed volumes, exported as one multi-solid STEP (plus three STLs for
rendering):

  winding   copper    -- volumetric heat source (copper loss)
  housing   aluminium -- finned sleeve, closed at both ends
  air       air       -- box with the motor subtracted

For CHT the fluid volume has to be built, not extracted from a solid part:
the air domain is the box minus the housing, so its walls are exactly the
wetted motor surfaces.

Why the housing is closed
-------------------------
Every material interface has exactly two materials on it:

    winding | housing     (cavity surface: lateral + both end discs)
    housing | air         (outer surface + fins)

An earlier version left the housing open at both ends with the winding
protruding through it. That put the winding's lateral surface in coincidence
with two different solids at once -- the housing bore along the sleeve and the
air along the overhang -- so winding, housing and air all met along a curve.
SimScale's contact detection does pairwise booleans and fails on that triple
line ("failure of boolean underlying the operation between A1_I7_B3_TE3,
A1_I11_B2_TE38"), which silently leaves the solid uncoupled. Closing the ends
removes the triple line entirely.

Consequence: the winding has no wetted area at all, so every watt leaves
through the housing. That is also what a potted servo stator does.

Axis convention: motor axis and flow direction are both +x.

Run:  python make_geometry.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import trimesh

# --------------------------------------------------------------------------
# Design parameters (metres)
# --------------------------------------------------------------------------
R_WIND = 26.0e-3         # winding radius = housing cavity radius
L_WIND = 112.0e-3        # winding length = finned sleeve length
R_SLEEVE_OUT = 32.0e-3   # housing outer radius under the fins
L_HOUSING = L_WIND       # finned sleeve length
N_FINS = 12
R_FIN_TIP = 42.0e-3      # fin tip radius -> 10 mm fin height
T_FIN = 1.5e-3           # fin thickness
T_END = 2.0e-3           # housing end wall thickness

# Fins must PENETRATE the sleeve, not sit tangent to it. A flat fin base at
# r = R_SLEEVE_OUT only touches the cylindrical sleeve along a line, so the
# boolean union leaves it as a separate body. Embedding them 2 mm into the
# sleeve fuses them.
EMBED = 2.0e-3

# Air domain
X_UP, X_DOWN = 150.0e-3, 450.0e-3   # upstream / downstream extent from centre
HALF_W = 150.0e-3                    # half width and half height

# Operating point
Q_WINDING = 80.0         # W, copper loss
T_INLET = 298.15         # K (25 C)
U_INLET = 3.0            # m/s

# Material properties at ~300 K
COPPER = {"k": 385.0, "rho": 8933.0, "cp": 385.0}
ALUMINIUM = {"k": 202.0, "rho": 2700.0, "cp": 900.0}
AIR = {"k": 0.0261, "rho": 1.18, "cp": 1005.0, "mu": 1.85e-5}
AIR["nu"] = AIR["mu"] / AIR["rho"]
AIR["pr"] = AIR["mu"] * AIR["cp"] / AIR["k"]

OUT = Path(__file__).resolve().parent


def cylinder_x(radius: float, length: float, centre_x: float = 0.0,
               sections: int = 128) -> trimesh.Trimesh:
    """Cylinder with its axis along +x."""
    m = trimesh.creation.cylinder(radius=radius, height=length, sections=sections)
    m.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2.0, [0, 1, 0]))
    m.apply_translation([centre_x, 0.0, 0.0])
    return m


def build_winding() -> trimesh.Trimesh:
    return cylinder_x(R_WIND, L_WIND)


def build_housing() -> trimesh.Trimesh:
    """Closed finned sleeve: outer cylinder, cavity, 12 axial fins, ONE solid."""
    outer = cylinder_x(R_SLEEVE_OUT, L_HOUSING + 2.0 * T_END)
    cavity = cylinder_x(R_WIND, L_HOUSING)
    parts = [trimesh.boolean.difference([outer, cavity])]

    # The fin root sits EMBED below the sleeve surface so the two solids
    # genuinely overlap: a fin base tangent to the sleeve only touches it along
    # a line and never fuses.
    r_root = R_SLEEVE_OUT - EMBED
    r_mid = 0.5 * (r_root + R_FIN_TIP)
    for i in range(N_FINS):
        fin = trimesh.creation.box(extents=[L_HOUSING, T_FIN, R_FIN_TIP - r_root])
        fin.apply_translation([0.0, 0.0, r_mid])
        fin.apply_transform(
            trimesh.transformations.rotation_matrix(2.0 * np.pi * i / N_FINS, [1, 0, 0])
        )
        parts.append(fin)

    return trimesh.boolean.union(parts)


def build_air(housing: trimesh.Trimesh, winding: trimesh.Trimesh) -> trimesh.Trimesh:
    """Fluid domain: the box minus the whole motor.

    Both solids are subtracted, and the winding matters even though it is
    enclosed. Subtracting only the housing leaves the housing cavity behind as
    a sealed air pocket -- a second, disconnected body inside the fluid domain.
    Worse, that pocket's wall is exactly the winding's surface, so air and
    winding would share a face and the triple line comes straight back. Cutting
    the winding out removes the pocket and leaves one connected air region
    whose only walls are the housing's outer surface and fins.
    """
    box = trimesh.creation.box(extents=[X_UP + X_DOWN, 2 * HALF_W, 2 * HALF_W])
    box.apply_translation([0.5 * (X_DOWN - X_UP), 0.0, 0.0])
    return trimesh.boolean.difference([box, housing, winding])


def build_step() -> dict:
    """Build the three solids as exact B-rep and export a multi-solid STEP.

    A single STL cannot carry this case. SimScale treats one imported file as
    one geometry, and a CHT simulation can only take a single geometry as its
    domain, so importing three STL files gave three geometries that could never
    be combined. A STEP file carries three independent solids, which is what
    the material assignment needs.

    Returns the measured volume of each solid, for checking against trimesh.
    """
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 1)
    occ = gmsh.model.occ

    try:
        # --- winding: cylinder along +x -----------------------------------
        winding = occ.addCylinder(-0.5 * L_WIND, 0, 0, L_WIND, 0, 0, R_WIND)

        # --- housing: closed sleeve + embedded fins ------------------------
        l_outer = L_HOUSING + 2.0 * T_END
        outer = occ.addCylinder(-0.5 * l_outer, 0, 0, l_outer, 0, 0, R_SLEEVE_OUT)
        cavity = occ.addCylinder(-0.5 * L_HOUSING, 0, 0, L_HOUSING, 0, 0, R_WIND)
        sleeve, _ = occ.cut([(3, outer)], [(3, cavity)])

        r_root = R_SLEEVE_OUT - EMBED
        fins = []
        for i in range(N_FINS):
            f = occ.addBox(-0.5 * L_HOUSING, -0.5 * T_FIN, r_root,
                           L_HOUSING, T_FIN, R_FIN_TIP - r_root)
            occ.rotate([(3, f)], 0, 0, 0, 1, 0, 0, 2.0 * np.pi * i / N_FINS)
            fins.append((3, f))
        housing, _ = occ.fuse(sleeve, fins)

        # --- air: box minus the whole motor ---------------------------------
        box = occ.addBox(-X_UP, -HALF_W, -HALF_W, X_UP + X_DOWN, 2 * HALF_W, 2 * HALF_W)
        # Copies are cut, not the originals: cut() deletes its tools by default
        # and the housing and winding still have to exist as their own solids.
        # The winding is cut as well so the housing cavity does not survive as
        # a sealed air pocket whose wall is the winding's surface.
        housing_tool = occ.copy(housing)
        winding_tool = occ.copy([(3, winding)])
        air, _ = occ.cut([(3, box)], housing_tool + winding_tool)

        occ.synchronize()

        tags = {"winding": winding, "housing": housing[0][1], "air": air[0][1]}
        if len(set(tags.values())) != 3:
            raise SystemExit(f"STEP build did not produce three distinct solids: {tags}")

        for name, tag in tags.items():
            gmsh.model.addPhysicalGroup(3, [tag], name=name)
        gmsh.model.occ.synchronize()

        # Mass with unit density is the volume, so this is a direct check that
        # each solid came out the right size rather than an empty shell.
        measured = {k: float(gmsh.model.occ.getMass(3, v)) for k, v in tags.items()}

        # The OCC writer always declares millimetres, and only rescales the
        # coordinates when it is told what unit the model is already in. With
        # the default ("MM") a model built in metres is written as 0.45 and
        # declared mm, so every consumer -- SimScale included -- sees a part
        # 1000x too small. Setting the source unit to M makes the writer scale
        # metres to millimetres, which matches the declaration.
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.write(str(OUT / "motor_assembly.step"))
    finally:
        gmsh.finalize()

    # Volumes only -- the unit check runs in main() so its result lands in the
    # report as a size, not as another volume that would be scaled to mL.
    verify_step_units(OUT / "motor_assembly.step", X_UP + X_DOWN)
    return measured


SI_PREFIX = {"": 1.0, "MILLI": 1e-3, "CENTI": 1e-2, "MICRO": 1e-6, "KILO": 1e3}


def verify_step_units(path: Path, expected_extent_m: float) -> float:
    """Read the STEP the way a consumer does and check the part is to scale.

    Nothing inside gmsh catches a unit mismatch: the writer declares one unit
    and writes coordinates in another, and a round-trip back through gmsh
    agrees with itself either way. So this parses the file directly -- declared
    length unit times written coordinates -- and refuses to ship a part whose
    real size is wrong. A 1000x-small domain still meshes and still converges;
    it just solves a different, physically meaningless case, which is exactly
    the kind of error that has to be caught before it reaches a report.
    """
    text = path.read_text(errors="ignore")

    unit = re.search(r"SI_UNIT\(\s*(?:\.([A-Z]+)\.|\$)\s*,\s*\.METRE\.\s*\)", text)
    if not unit:
        raise SystemExit("no length unit found in the STEP header")
    prefix = unit.group(1) or ""
    if prefix not in SI_PREFIX:
        raise SystemExit(f"unhandled STEP length prefix {prefix!r}")
    scale = SI_PREFIX[prefix]

    # Only 3D points are model space. A STEP file also carries 2D points in
    # (u, v) parameter space for the pcurves of every edge, and those routinely
    # sit outside the part -- picking one up makes the extent read as nonsense
    # (the first run of this check reported 750 mm because of a pcurve point at
    # u = 600), so the regex insists on all three coordinates.
    xs = [
        float(m.group(1))
        for m in re.finditer(
            r"CARTESIAN_POINT\s*\(\s*'[^']*'\s*,\s*\(\s*([-0-9.Ee+]+)\s*,"
            r"\s*([-0-9.Ee+]+)\s*,\s*([-0-9.Ee+]+)\s*\)",
            text,
        )
    ]
    if not xs:
        raise SystemExit("no 3D CARTESIAN_POINT coordinates found in the STEP")
    extent = (max(xs) - min(xs)) * scale
    if abs(extent - expected_extent_m) > 0.01 * expected_extent_m:
        raise SystemExit(
            f"STEP is not to scale: x reads {extent * 1e3:.1f} mm but should be "
            f"{expected_extent_m * 1e3:.1f} mm. Declared unit is "
            f"{prefix or 'METRE'} and the coordinates run {min(xs):.4g} to "
            f"{max(xs):.4g}, so a consumer scales the part by {scale:g}."
        )
    return extent


def analytics(housing: trimesh.Trimesh, winding: trimesh.Trimesh) -> dict:
    # --- interface and wetted areas -------------------------------------
    # The winding is enclosed, so it contacts the housing over its whole
    # surface (lateral + both end discs) and has no wetted area at all.
    a_interface = 2.0 * np.pi * R_WIND * L_HOUSING + 2.0 * np.pi * R_WIND ** 2
    a_housing_ext = float(housing.area) - a_interface
    a_winding_ext = 0.0
    a_total = a_housing_ext

    # --- air-side heat transfer coefficient -----------------------------
    # Axial flow along the fins. A flat-plate estimate over the motor length
    # is the honest anchor here: Churchill-Bernstein is a crossflow
    # correlation and does not apply to flow along the axis.
    re_l = U_INLET * L_WIND / AIR["nu"]
    nu_l = 0.037 * re_l ** 0.8 * AIR["pr"] ** (1.0 / 3.0)   # turbulent flat plate
    h = nu_l * AIR["k"] / L_WIND

    # --- fin efficiency --------------------------------------------------
    m = np.sqrt(2.0 * h / (ALUMINIUM["k"] * T_FIN))
    ml = m * (R_FIN_TIP - R_SLEEVE_OUT)
    eta_fin = np.tanh(ml) / ml

    a_fins = N_FINS * (2.0 * (R_FIN_TIP - R_SLEEVE_OUT) * L_HOUSING + T_FIN * L_HOUSING)
    a_base = a_total - a_fins
    a_eff = a_base + eta_fin * a_fins
    r_conv = 1.0 / (h * a_eff)

    # --- conduction through the housing wall -----------------------------
    r_wall = np.log(R_SLEEVE_OUT / R_WIND) / (2.0 * np.pi * ALUMINIUM["k"] * L_HOUSING)

    # --- lumped prediction ----------------------------------------------
    dt_conv = Q_WINDING * r_conv
    dt_air_bulk = Q_WINDING / (AIR["rho"] * U_INLET * (2 * HALF_W) ** 2 * AIR["cp"])
    q_vol = Q_WINDING / float(winding.volume)
    # uniform generation in a copper cylinder: centre-to-surface drop
    dt_winding_internal = q_vol * R_WIND ** 2 / (4.0 * COPPER["k"])

    return {
        "geometry": {
            "winding_radius_m": R_WIND,
            "winding_length_m": L_WIND,
            "housing_outer_radius_m": R_SLEEVE_OUT,
            "fin_tip_radius_m": R_FIN_TIP,
            "fin_height_m": R_FIN_TIP - R_SLEEVE_OUT,
            "fin_thickness_m": T_FIN,
            "n_fins": N_FINS,
            "housing_length_m": L_HOUSING + 2.0 * T_END,
            "end_wall_thickness_m": T_END,
            "interface_area_m2": a_interface,
            "housing_wetted_area_m2": a_housing_ext,
            "winding_exposed_area_m2": a_winding_ext,
            "total_wetted_area_m2": a_total,
            "fin_area_m2": a_fins,
            "air_domain_m3": (X_UP + X_DOWN) * (2 * HALF_W) ** 2,
        },
        "operating_point": {
            "heat_load_W": Q_WINDING,
            "heat_source_W_m3": q_vol,
            "inlet_temperature_K": T_INLET,
            "inlet_velocity_m_s": U_INLET,
            "copper": COPPER,
            "aluminium": ALUMINIUM,
            "air": AIR,
        },
        "prediction": {
            "re_length": re_l,
            "nusselt_flat_plate": nu_l,
            "h_W_m2K": h,
            "fin_parameter_m": m,
            "fin_efficiency": eta_fin,
            "effective_area_m2": a_eff,
            "convective_resistance_K_W": r_conv,
            "wall_conduction_resistance_K_W": r_wall,
            "housing_over_air_K": dt_conv,
            "predicted_housing_temperature_K": T_INLET + dt_conv,
            "predicted_bore_temperature_K": T_INLET + Q_WINDING * (r_conv + r_wall),
            "winding_internal_drop_K": dt_winding_internal,
            "air_bulk_rise_K": dt_air_bulk,
            "note": (
                "h from a turbulent flat-plate correlation over the motor length, "
                "not Churchill-Bernstein: the flow is axial, not crossflow. "
                "The lumped resistance assumes the housing is isothermal, which "
                "it is not; it is a scale check, not a prediction of T_max. "
                "The winding is enclosed, so all 80 W must cross the air-side "
                "resistance; wall conduction is included as a separate term."
            ),
        },
    }


def main() -> None:
    winding = build_winding()
    housing = build_housing()
    air = build_air(housing, winding)

    for name, mesh in [("winding", winding), ("housing", housing), ("air", air)]:
        if not mesh.is_watertight:
            raise SystemExit(f"{name} is not watertight -- refusing to export")
        # A part may legitimately have several disconnected SURFACES: a solid
        # with an internal cavity has an outer skin plus the cavity, and the
        # cavity carries negative volume. What must be unique is the number of
        # material bodies, i.e. positive-volume shells. SimScale counts each
        # loose solid as a separate volume, so two of them silently turns a
        # part into an unusable multi-region geometry.
        shells = mesh.split(only_watertight=False)
        solids = [s for s in shells if s.volume > 0]
        if len(solids) != 1:
            raise SystemExit(
                f"{name} has {len(solids)} material bodies -- each part must be a "
                f"single solid, or the CHT region assignment breaks"
            )

    winding.export(OUT / "motor_winding.stl")
    housing.export(OUT / "motor_housing.stl")
    air.export(OUT / "air_domain.stl")

    # The STLs above are for rendering. SimScale needs the three-solid STEP;
    # see build_step() for why.
    step_volumes = build_step()

    report = {
        "parts": {
            "winding": {
                "material": "copper",
                "volume_mL": float(winding.volume) * 1e6,
                "area_m2": float(winding.area),
                "bounds_m": winding.bounds.tolist(),
                "n_faces": int(len(winding.faces)),
            },
            "housing": {
                "material": "aluminium",
                "volume_mL": float(housing.volume) * 1e6,
                "area_m2": float(housing.area),
                "bounds_m": housing.bounds.tolist(),
                "n_faces": int(len(housing.faces)),
            },
            "air": {
                "material": "air",
                "volume_mL": float(air.volume) * 1e6,
                "area_m2": float(air.area),
                "bounds_m": air.bounds.tolist(),
                "n_faces": int(len(air.faces)),
            },
            "motor": {
                "volume_mL": (float(housing.volume) + float(winding.volume)) * 1e6,
                "note": "housing + winding, reported separately; not unioned",
            },
            "assembly": {
                "file": "motor_assembly.step",
                "n_solids": len(step_volumes),
                "x_extent_m": X_UP + X_DOWN,
                "volumes_mL": {k: round(v * 1e6, 2) for k, v in step_volumes.items()},
                "note": (
                    "three independent B-rep solids in one STEP, because a CHT "
                    "simulation takes a single geometry and the winding/housing "
                    "interface is a shared face that an STL cannot keep separate"
                ),
            },
        },
        **analytics(housing, winding),
    }
    (OUT / "design_point.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8", newline="\n"
    )

    p = report["parts"]
    print(f"winding  {p['winding']['volume_mL']:8.3f} mL   {p['winding']['area_m2']*1e4:7.2f} cm2")
    print(f"housing  {p['housing']['volume_mL']:8.3f} mL   {p['housing']['area_m2']*1e4:7.2f} cm2")
    print(f"air      {p['air']['volume_mL']:8.1f} mL   {p['air']['area_m2']*1e4:7.1f} cm2")
    print()
    pr = report["prediction"]
    print(f"Re_L              {pr['re_length']:10.0f}")
    print(f"h                 {pr['h_W_m2K']:10.1f} W/m2K")
    print(f"fin efficiency    {pr['fin_efficiency']:10.3f}")
    print(f"effective area    {pr['effective_area_m2']*1e4:10.2f} cm2")
    print(f"R_conv            {pr['convective_resistance_K_W']:10.3f} K/W")
    print(f"R_wall            {pr['wall_conduction_resistance_K_W']:10.4f} K/W")
    print(f"housing over air  {pr['housing_over_air_K']:10.1f} K")
    print(f"predicted housing {pr['predicted_housing_temperature_K'] - 273.15:10.1f} C")
    print(f"predicted bore    {pr['predicted_bore_temperature_K'] - 273.15:10.1f} C")
    print(f"winding internal  {pr['winding_internal_drop_K']:10.3f} K")
    print(f"air bulk rise     {pr['air_bulk_rise_K']:10.3f} K")
    print(f"\nheat source       {report['operating_point']['heat_source_W_m3']:.3e} W/m3")
    print(f"Wrote 3 STLs + STEP + design_point.json to {OUT}")


if __name__ == "__main__":
    main()
