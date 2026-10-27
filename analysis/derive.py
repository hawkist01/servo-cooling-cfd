"""
Derive validation metrics from the SimScale conjugate-heat-transfer run.

Why this file exists
--------------------
Every measured number is read from ``results/measurements.json``, which records
what the solver produced and where it came from.  Every reference value comes
from ``geometry/design_point.json``, written before the mesh was generated.
Nothing here invents a number.

Where the solid temperatures come from
--------------------------------------
The two solid surface-data result controls produced no plots in this solver
version, so there are no surface-averaged solid values.  The solid temperatures
therefore come from the result archive itself: ``result.zip`` holds OpenFOAM
ASCII fields per region, and ``1000/region2/T.gz`` / ``1000/region3/T.gz`` were
decompressed and parsed.  Those are cell-centre values, so every solid average
below is a **cell (volume) mean**, not an area mean.  That distinction is stated
in the report; it costs about 0.1 K here because the whole housing spans 2.8 K.

The wetted-surface temperature
------------------------------
The winding is fully enclosed, so the housing carries both the wetted faces and
the bore patch that takes heat from the winding.  The two are separated by known
geometry:

    A_wetted = A_housing - A_bore

    T_wetted = (A_housing * T_housing - A_bore * T_winding) / A_wetted

Both areas come from the CAD.  The correction moves the answer by ~0.4 K, which
is itself the point: the solid is nearly isothermal.

The two energy balances
-----------------------
The area-average outlet temperature is **not** the mixed-mean (bulk) outlet
temperature: it weights every unit of area equally, so the slow warm wake counts
as much as the fast cold core.  Both are reported.  The gap between them is the
velocity-temperature covariance, which only field data could close.

Usage
-----
    python derive.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
MEASUREMENTS = RESULTS / "measurements.json"
DESIGN_POINT = ROOT / "geometry" / "design_point.json"

Q_IMPOSED_W = 80.0
G = 9.81
T_IN_NOMINAL = 298.15


def load() -> tuple[dict, dict]:
    if not MEASUREMENTS.exists():
        raise SystemExit(f"missing {MEASUREMENTS}")
    meas = json.loads(MEASUREMENTS.read_text(encoding="utf-8"))
    dp = json.loads(DESIGN_POINT.read_text(encoding="utf-8"))
    return meas, dp


def main() -> None:
    meas, dp = load()

    pred = dp["prediction"]
    geom = dp["geometry"]
    op = dp["operating_point"]

    # --- fluid properties actually used by the solver ---------------------
    air = meas["fluid_properties"]
    rho = air["rho_kg_m3"]
    cp = air["cp_J_kgK"]
    k_air = air["k_W_mK"]
    nu = air["nu_m2_s"]
    pr = air["prandtl"]

    solids = meas["solid_properties"]
    k_al = solids["aluminium"]["k_W_mK"]

    bc = meas["boundary_conditions"]
    u_in = bc["inlet_velocity_m_s"]
    a_in = bc["inlet_area_m2"]

    rc = meas["result_controls"]
    t_in = rc["inlet_T_K"]
    t_out = rc["outlet_T_K"]

    # --- mass flow and the two bulk rises ---------------------------------
    mdot = rho * u_in * a_in
    mcp = mdot * cp

    dt_area = t_out - t_in
    q_area = mcp * dt_area

    dt_mixed = Q_IMPOSED_W / mcp
    area_bias_pct = 100.0 * (dt_area - dt_mixed) / dt_mixed

    mass_err_pct = 100.0 * (rc["outlet_Ux_m_s"] / rc["inlet_Ux_m_s"] - 1.0)

    # --- solid temperatures from the parsed fields ------------------------
    sf = meas["solid_fields"]
    housing = sf["housing"]
    winding = sf["winding"]

    t_housing = housing["mean_K"]
    t_winding = winding["mean_K"]
    t_winding_max = winding["max_K"]
    housing_spread = housing["max_K"] - housing["min_K"]
    winding_spread = winding["max_K"] - winding["min_K"]

    a_bore = geom["interface_area_m2"]
    a_wetted = geom["housing_wetted_area_m2"]
    a_housing = a_wetted + a_bore

    t_wetted = (a_housing * t_housing - a_bore * t_winding) / a_wetted

    # --- thermal path ------------------------------------------------------
    t_bulk_mean = 0.5 * (t_in + t_out)
    superheat = t_wetted - t_bulk_mean
    metal_drop = t_winding - t_wetted

    # --- effective heat transfer coefficient ------------------------------
    a_eff = pred["effective_area_m2"]
    h_eff = Q_IMPOSED_W / (a_eff * superheat)
    h_pred = pred["h_W_m2K"]
    h_err_pct = 100.0 * (h_eff - h_pred) / h_pred

    # Same correlation re-evaluated with the solver's own air properties and
    # the wetted length, so the property and length choices are visible rather
    # than buried in the comparison.
    length_wetted = geom["housing_length_m"]
    re_l = u_in * length_wetted / nu
    nu_eff = h_eff * length_wetted / k_air
    nu_pred_at_l = 0.037 * re_l**0.8 * pr ** (1.0 / 3.0)
    h_pred_same_props = nu_pred_at_l * k_air / length_wetted

    # --- buoyancy ratio ----------------------------------------------------
    beta = air.get("beta_1_K", 1.0 / t_in)
    gr_l = G * beta * (t_wetted - t_in) * length_wetted**3 / nu**2
    ri = gr_l / re_l**2

    # --- thermal resistance ------------------------------------------------
    r_th_cfd = (t_wetted - t_in) / Q_IMPOSED_W
    r_th_pred = pred["convective_resistance_K_W"] + pred["wall_conduction_resistance_K_W"]

    # --- Biot numbers ------------------------------------------------------
    bi_fin = h_eff * geom["fin_height_m"] / k_al
    bi_wall = h_eff * geom["end_wall_thickness_m"] / k_al

    # --- what radiation would add -----------------------------------------
    # Not modelled (enableRadiation is false).  Quantify the omission so the
    # reader can judge it: the black-body flux to a 25 C surround at the
    # measured surface temperature, as a fraction of the convective flux, for
    # bare machined aluminium (eps ~ 0.1) and for anodised/painted metal
    # (eps ~ 0.8).
    sigma = 5.670374419e-8
    q_rad_bb = sigma * (t_wetted**4 - t_in**4)
    q_conv_flux = Q_IMPOSED_W / a_eff
    rad_frac_bare = 100.0 * 0.1 * q_rad_bb / q_conv_flux
    rad_frac_coated = 100.0 * 0.8 * q_rad_bb / q_conv_flux

    # --- prediction comparators -------------------------------------------
    housing_over_inlet = t_wetted - t_in
    housing_over_inlet_pred = pred["housing_over_air_K"]

    hot_spot_pred = pred["predicted_bore_temperature_K"] + pred["winding_internal_drop_K"]
    hot_spot_over_inlet = t_winding_max - t_in
    hot_spot_pred_over_inlet = hot_spot_pred - T_IN_NOMINAL

    rows = {
        "mdot_kg_s": mdot,
        "mcp_W_K": mcp,
        "dt_area_K": dt_area,
        "q_area_W": q_area,
        "dt_mixed_K": dt_mixed,
        "area_vs_mixed_bias_pct": area_bias_pct,
        "mass_conservation_err_pct": mass_err_pct,
        "t_in_K": t_in,
        "t_out_K": t_out,
        "t_bulk_mean_K": t_bulk_mean,
        "t_housing_avg_K": t_housing,
        "t_winding_avg_K": t_winding,
        "t_wetted_avg_K": t_wetted,
        "t_winding_max_K": t_winding_max,
        "t_housing_min_K": housing["min_K"],
        "t_housing_max_K": housing["max_K"],
        "housing_spread_K": housing_spread,
        "winding_spread_K": winding_spread,
        "winding_peak_over_mean_K": t_winding_max - t_winding,
        "metal_drop_K": metal_drop,
        "superheat_K": superheat,
        "housing_over_inlet_K": housing_over_inlet,
        "housing_over_inlet_predicted_K": housing_over_inlet_pred,
        "housing_vs_prediction_pct": 100.0
        * (housing_over_inlet - housing_over_inlet_pred)
        / housing_over_inlet_pred,
        "hot_spot_over_inlet_K": hot_spot_over_inlet,
        "hot_spot_predicted_over_inlet_K": hot_spot_pred_over_inlet,
        "hot_spot_vs_prediction_pct": 100.0
        * (hot_spot_over_inlet - hot_spot_pred_over_inlet)
        / hot_spot_pred_over_inlet,
        "h_effective_W_m2K": h_eff,
        "h_predicted_W_m2K": h_pred,
        "h_predicted_same_props_W_m2K": h_pred_same_props,
        "h_vs_prediction_pct": h_err_pct,
        "h_vs_same_props_pct": 100.0 * (h_eff - h_pred_same_props) / h_pred_same_props,
        "nusselt_effective": nu_eff,
        "nusselt_predicted": pred["nusselt_flat_plate"],
        "nusselt_predicted_at_L": nu_pred_at_l,
        "reynolds_length": re_l,
        "reynolds_length_predicted": pred["re_length"],
        "grashof_length": gr_l,
        "richardson_number": ri,
        "richardson_percent": 100.0 * ri,
        "prandtl": pr,
        "thermal_resistance_cfd_K_W": r_th_cfd,
        "thermal_resistance_predicted_K_W": r_th_pred,
        "thermal_resistance_vs_prediction_pct": 100.0 * (r_th_cfd - r_th_pred) / r_th_pred,
        "biot_fin": bi_fin,
        "biot_wall": bi_wall,
        "fin_efficiency": pred["fin_efficiency"],
        "radiative_flux_blackbody_W_m2": q_rad_bb,
        "convective_flux_W_m2": q_conv_flux,
        "radiation_share_bare_pct": rad_frac_bare,
        "radiation_share_coated_pct": rad_frac_coated,
        "wetted_area_m2": a_wetted,
        "effective_area_m2": a_eff,
        "length_wetted_m": length_wetted,
    }

    # ---- console report --------------------------------------------------
    line = "=" * 78
    print(line)
    print("SERVO MOTOR - CONJUGATE HEAT TRANSFER")
    print(line)
    print(f"mesh                 {meas['mesh']['cells']:,} cells / "
          f"{meas['mesh']['nodes']:,} nodes")
    print(f"air                  rho {rho} kg/m3, cp {cp} J/kgK, "
          f"nu {nu:.4e} m2/s, Pr {pr}")
    print()
    print("ENERGY BALANCE")
    print(f"   mdot              {mdot:.5f} kg/s   (rho U A)")
    print(f"   mdot cp           {mcp:.3f} W/K")
    print(f"   T_in              {t_in:9.4f} K")
    print(f"   T_out (area avg)  {t_out:9.4f} K    dT = {dt_area:.4f} K")
    print(f"   -> area-average Q {q_area:9.3f} W    ({area_bias_pct:+.1f} % vs 80 W)")
    print(f"   mixed-mean dT     {dt_mixed:9.4f} K    -> Q = {Q_IMPOSED_W:.1f} W")
    print(f"   mass conservation {mass_err_pct:+.2e} %  (outlet/inlet Ux)")
    print()
    print("TEMPERATURES (cell means of the parsed solid fields)")
    print(f"   housing (all)     {t_housing:9.4f} K  ({t_housing - 273.15:7.3f} C)"
          f"   span {housing_spread:.3f} K")
    print(f"   winding           {t_winding:9.4f} K  ({t_winding - 273.15:7.3f} C)"
          f"   span {winding_spread:.3f} K")
    print(f"   winding peak      {t_winding_max:9.4f} K  ({t_winding_max - 273.15:7.3f} C)"
          f"   <- measured, not estimated")
    print(f"   wetted surface    {t_wetted:9.4f} K  ({t_wetted - 273.15:7.3f} C)")
    print()
    print("THERMAL PATH")
    print(f"   air bulk rise     {dt_mixed:9.4f} K")
    print(f"   air film          {superheat:9.4f} K   <- the bottleneck")
    print(f"   metal             {metal_drop:9.4f} K   (wetted surface -> winding)")
    print()
    print("HEAT TRANSFER COEFFICIENT")
    print(f"   h effective       {h_eff:9.3f} W/m2K  vs {h_pred:.3f} predicted  "
          f"({h_err_pct:+.1f} %)")
    print(f"   h (same props/L)  {h_pred_same_props:9.3f} W/m2K  <- correlation on the "
          f"solver's air")
    print(f"   Nu at L={length_wetted*1000:.0f} mm {nu_eff:9.2f}       vs "
          f"{nu_pred_at_l:.2f} (correlation, same L)")
    print(f"   Re_L              {re_l:9.0f}       (prediction used "
          f"{pred['re_length']:.0f} at L={geom['winding_length_m']*1000:.0f} mm)")
    print(f"   Gr_L              {gr_l:9.3e}   Ri = {ri:.4f}")
    print()
    print("THERMAL RESISTANCE")
    print(f"   R = (T_w - T_in)/Q {r_th_cfd:8.4f} K/W  vs {r_th_pred:.4f}  "
          f"({rows['thermal_resistance_vs_prediction_pct']:+.1f} %)")
    print(f"   Bi (fin)          {bi_fin:9.2e}   Bi (end wall) {bi_wall:.2e}")
    print()
    print("RADIATION (not modelled)")
    print(f"   black-body flux   {q_rad_bb:9.1f} W/m2  vs {q_conv_flux:.0f} W/m2 convective")
    print(f"   share, eps=0.1    {rad_frac_bare:9.1f} %     eps=0.8  {rad_frac_coated:.1f} %")

    out = RESULTS / "summary.json"
    out.write_text(
        json.dumps(
            {
                "project": meas["project"],
                "mesh": meas["mesh"],
                "fluid_properties": air,
                "solid_properties": solids,
                "boundary_conditions": bc,
                "result_controls": rc,
                "solid_fields": sf,
                "convergence": meas.get("convergence"),
                "measured": rows,
                "prediction": {
                    "h_W_m2K": h_pred,
                    "nusselt_flat_plate": pred["nusselt_flat_plate"],
                    "re_length": pred["re_length"],
                    "effective_area_m2": a_eff,
                    "housing_over_air_K": housing_over_inlet_pred,
                    "air_bulk_rise_K": pred["air_bulk_rise_K"],
                    "fin_efficiency": pred["fin_efficiency"],
                    "fin_parameter_m": pred["fin_parameter_m"],
                    "thermal_resistance_K_W": r_th_pred,
                    "convective_resistance_K_W": pred["convective_resistance_K_W"],
                    "wall_conduction_resistance_K_W": pred[
                        "wall_conduction_resistance_K_W"
                    ],
                    "predicted_housing_temperature_K": pred[
                        "predicted_housing_temperature_K"
                    ],
                    "predicted_bore_temperature_K": pred["predicted_bore_temperature_K"],
                    "winding_internal_drop_K": pred["winding_internal_drop_K"],
                    "predicted_hot_spot_K": hot_spot_pred,
                    "wetted_area_m2": geom["total_wetted_area_m2"],
                    "fin_area_m2": geom["fin_area_m2"],
                    "interface_area_m2": geom["interface_area_m2"],
                    "housing_area_m2": a_housing,
                    "n_fins": geom["n_fins"],
                    "fin_height_m": geom["fin_height_m"],
                    "fin_thickness_m": geom["fin_thickness_m"],
                    "housing_length_m": geom["housing_length_m"],
                    "winding_length_m": geom["winding_length_m"],
                    "winding_radius_m": geom["winding_radius_m"],
                    "housing_outer_radius_m": geom["housing_outer_radius_m"],
                    "fin_tip_radius_m": geom["fin_tip_radius_m"],
                    "end_wall_thickness_m": geom["end_wall_thickness_m"],
                    "air_domain_m3": geom["air_domain_m3"],
                    "heat_load_W": op["heat_load_W"],
                    "heat_source_W_m3": op["heat_source_W_m3"],
                    "inlet_temperature_K": op["inlet_temperature_K"],
                    "inlet_velocity_m_s": op["inlet_velocity_m_s"],
                    "air_properties_used": op["air"],
                },
                "provenance": meas.get("provenance"),
            },
            indent=2,
        ),
        encoding="utf-8",
        newline="\n",
    )
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
