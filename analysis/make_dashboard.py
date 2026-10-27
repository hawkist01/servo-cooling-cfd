"""
Build docs/dashboard.html - a single self-contained summary page.

Reads results/summary.json and inlines the figures as base64 so the file opens
offline with no external assets.

Run:  python make_dashboard.py
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
FIG = ROOT / "figures"
DOCS = ROOT / "docs"

CSS = """
:root {
  --ink: #12212e; --muted: #5b6b7a; --line: #d9e0e6;
  --blue: #1f4e79; --orange: #c55a11; --green: #2e7d32; --copper: #b87333;
  --bg: #f6f8fa; --card: #ffffff;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
       font: 15px/1.55 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
.wrap { max-width: 1080px; margin: 0 auto; padding: 40px 24px 72px; }
header { border-bottom: 3px solid var(--blue); padding-bottom: 18px; margin-bottom: 30px; }
h1 { font-size: 26px; margin: 0 0 6px; letter-spacing: -0.3px; }
.sub { color: var(--muted); font-size: 14px; }
h2 { font-size: 15px; text-transform: uppercase; letter-spacing: 1.1px;
     color: var(--muted); margin: 38px 0 14px; font-weight: 600; }
.kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; }
.kpi { background: var(--card); border: 1px solid var(--line); border-radius: 10px;
       padding: 16px 18px; }
.kpi .v { font-size: 27px; font-weight: 650; letter-spacing: -0.5px; }
.kpi .l { font-size: 12px; color: var(--muted); margin-top: 3px; line-height: 1.35; }
.kpi .d { font-size: 12px; margin-top: 6px; color: var(--muted); }
.blue { color: var(--blue); } .orange { color: var(--orange); }
.green { color: var(--green); } .copper { color: var(--copper); }
figure { margin: 0 0 22px; background: var(--card); border: 1px solid var(--line);
         border-radius: 10px; padding: 14px; }
figure img { width: 100%; display: block; }
figcaption { font-size: 12.5px; color: var(--muted); margin-top: 10px; }
table { width: 100%; border-collapse: collapse; background: var(--card);
        border: 1px solid var(--line); border-radius: 10px; overflow: hidden;
        font-size: 14px; }
th, td { text-align: left; padding: 9px 14px; border-bottom: 1px solid var(--line); }
th { background: #eef2f5; font-size: 12px; text-transform: uppercase;
     letter-spacing: 0.6px; color: var(--muted); font-weight: 600; }
tr:last-child td { border-bottom: none; }
td.num { font-variant-numeric: tabular-nums; text-align: right; }
.note { background: #fff8e6; border: 1px solid #f0dca8; border-left: 4px solid var(--orange);
        border-radius: 8px; padding: 14px 18px; font-size: 14px; }
.note b { color: var(--orange); }
.grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
@media (max-width: 760px) { .grid2 { grid-template-columns: 1fr; } }
footer { margin-top: 44px; padding-top: 16px; border-top: 1px solid var(--line);
         font-size: 12.5px; color: var(--muted); }
a { color: var(--blue); }
"""


def img(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def main() -> None:
    summary_path = RESULTS / "summary.json"
    if not summary_path.exists():
        raise SystemExit(f"missing {summary_path} - run derive.py first")
    d = json.loads(summary_path.read_text(encoding="utf-8"))

    m = d["measured"]
    p = d["prediction"]
    mesh = d["mesh"]
    proj = d["project"]

    def kpi(value: str, label: str, colour: str = "", delta: str = "") -> str:
        cls = f' class="{colour}"' if colour else ""
        d_html = f'<div class="d">{delta}</div>' if delta else ""
        return (f'<div class="kpi"><div class="v"{cls}>{value}</div>'
                f'<div class="l">{label}</div>{d_html}</div>')

    kpis = "".join([
        kpi(f"{m['t_winding_max_K'] - 273.15:.1f} °C", "winding peak (measured)",
            "copper", f"predicted {p['predicted_hot_spot_K'] - 273.15:.1f} °C"),
        kpi(f"{m['housing_over_inlet_K']:.1f} K", "wetted surface over inlet air",
            "orange", f"{m['housing_vs_prediction_pct']:+.1f} % vs prediction"),
        kpi(f"{m['h_effective_W_m2K']:.1f}", "effective h [W/m²K]",
            "blue", f"predicted {p['h_W_m2K']:.1f}"),
        kpi(f"{m['superheat_K']:.1f} K", "of that, the air film",
            "", f"metal only {m['metal_drop_K']:.2f} K"),
        kpi(f"{m['area_vs_mixed_bias_pct']:+.1f} %", "area vs mixed-mean outlet T",
            "green", f"not an energy imbalance"),
    ])

    rows = [
        ("Wetted area", f"{m['wetted_area_m2']*1e4:.1f} cm²", "fins + barrel + end caps"),
        ("Fin efficiency", f"{p['fin_efficiency']:.4f}", "fins are isothermal"),
        ("Effective area", f"{m['effective_area_m2']*1e4:.1f} cm²", "η-weighted"),
        ("Re_L", f"{m['reynolds_length']:,.0f}", f"at L = {m['length_wetted_m']} m"),
        ("Nu, correlation", f"{m['nusselt_predicted_at_L']:.1f}", "0.037 Re^0.8 Pr^1/3, same L"),
        ("Nu, CFD", f"{m['nusselt_effective']:.1f}", "from measured superheat"),
        ("Convective R", f"{p['convective_resistance_K_W']:.4f} K/W", "99.8 % of the total"),
        ("Wall conduction R", f"{p['wall_conduction_resistance_K_W']:.5f} K/W", "0.2 %"),
        ("Biot (fin)", f"{m['biot_fin']:.2e}", "solid is lumped"),
        ("Ri = Gr/Re²", f"{m['richardson_number']:.4f}", "buoyancy omitted, justified"),
        ("Air bulk rise", f"{m['dt_mixed_K']:.4f} K", f"predicted {p['air_bulk_rise_K']:.4f} K"),
        ("Mass conservation", f"{m['mass_conservation_err_pct']:.2e} %", "outlet / inlet Ux"),
    ]
    table = "".join(
        f"<tr><td>{a}</td><td class='num'>{b}</td><td>{c}</td></tr>"
        for a, b, c in rows
    )

    conv_fig = ""
    if (FIG / "convergence.png").exists():
        conv_fig = f"""
<figure>
  <img alt="Convergence of monitors and residuals" src="{img(FIG / 'convergence.png')}">
  <figcaption>
    {d['convergence']['iterations'] if d.get('convergence') else '—'} iterations to steady
    state; the surface monitors stop moving long before the residual floor is reached.
  </figcaption>
</figure>"""

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Servo Motor CHT — CFD Summary</title>
<style>{CSS}</style>
</head>
<body>
<div class="wrap">

<header>
  <h1>Servo Motor — Conjugate Heat Transfer</h1>
  <div class="sub">
    80 W copper loss · 12-fin aluminium housing · 3 m/s axial air at 25 °C ·
    {mesh['cells']:,}-cell mesh · k-ω SST steady RANS
  </div>
</header>

<div class="kpis">{kpis}</div>

<h2>The thermal path</h2>
<figure>
  <img alt="Thermal path from air inlet to winding mean"
       src="{img(FIG / 'thermal_path.png')}">
  <figcaption>
    Every kelvin between the inlet air and the winding mean, split by the
    resistance that produces it. The metal accounts for
    {m['metal_drop_K']:.2f} K; the air film accounts for {m['superheat_K']:.1f} K.
  </figcaption>
</figure>

<h2>Prediction vs CFD</h2>
<div class="grid2">
  <figure>
    <img alt="Prediction vs CFD" src="{img(FIG / 'validation_prediction.png')}">
    <figcaption>Pre-run lumped model against the coupled solution.</figcaption>
  </figure>
  <figure>
    <img alt="Energy balance" src="{img(FIG / 'energy_balance.png')}">
    <figcaption>
      The area-average outlet temperature implies
      {m['q_area_W']:.1f} W of the 80 W imposed. The gap is the
      area-versus-mass weighting of the wake, not a conservation error.
    </figcaption>
  </figure>
</div>

{conv_fig}

<h2>Numbers behind it</h2>
<table>
  <tr><th>Quantity</th><th style="text-align:right">Value</th><th>Note</th></tr>
  {table}
</table>

<h2>The finding</h2>
<div class="note">
  <b>The winding has no wetted area.</b> The housing is closed at both ends, so
  all 80 W must cross the air-side film on the housing and fins. With
  Bi ≈ {m['biot_fin']:.0e} the solid is lumped, and the fin efficiency is
  {p['fin_efficiency']:.3f}, so the winding, bore and fin tips all sit inside a
  single {m['housing_spread_K']:.1f} K band. Only two levers move the answer:
  wetted area, or air velocity. Thicker walls and better copper do nothing.
</div>

<footer>
  SimScale project <a href="{proj['simscale_project_url']}">{proj['simscale_project_url']}</a><br>
  {proj['analysis_type']} · {proj['turbulence']} · {proj['time_dependency']}<br>
  Prediction fixed before the mesh was generated. Single mesh — no mesh
  independence study. Perfect thermal contact at all 66 interfaces, so the
  winding temperature is a lower bound.
</footer>

</div>
</body>
</html>
"""
    DOCS.mkdir(exist_ok=True)
    out = DOCS / "dashboard.html"
    out.write_text(html, encoding="utf-8", newline="\n")
    print(f"wrote {out}  ({len(html):,} chars)")


if __name__ == "__main__":
    main()
