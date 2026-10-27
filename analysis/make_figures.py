"""
Figures for the servo-motor conjugate-heat-transfer report.

Reads results/summary.json (written by derive.py) and, if present,
results/convergence.json, and produces:

  figures/thermal_path.png           where the 63 K of temperature rise comes from
  figures/validation_prediction.png  pre-run prediction vs CFD
  figures/energy_balance.png         area-average vs mixed-mean outlet temperature
  figures/convergence.png            monitors and residuals against iteration

Run:  python make_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
FIG = ROOT / "figures"

BLUE = "#1f4e79"
ORANGE = "#c55a11"
GREY = "#7f7f7f"
GREEN = "#2e7d32"
RED = "#b3261e"
COPPER = "#b87333"


def load() -> dict:
    p = RESULTS / "summary.json"
    if not p.exists():
        raise SystemExit(f"missing {p} - run derive.py first")
    return json.loads(p.read_text(encoding="utf-8"))


def style(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.set_axisbelow(True)


def fig_thermal_path(d: dict) -> None:
    """The ladder, and the share each step owns."""
    m = d["measured"]

    stages = [
        ("air bulk rise", m["dt_mixed_K"], BLUE),
        ("air film\n(wetted surface - bulk)", m["superheat_K"], ORANGE),
        ("metal\n(winding - wetted surface)", m["metal_drop_K"], COPPER),
    ]
    total = sum(s[1] for s in stages)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.5), dpi=160,
                             gridspec_kw={"width_ratios": [1.05, 1.0]})

    # --- left: share of the total rise -----------------------------------
    ax = axes[0]
    left = 0.0
    for label, value, colour in stages:
        share = 100.0 * value / total
        ax.barh(0, share, left=left, color=colour, height=0.42,
                edgecolor="white", linewidth=1.2, zorder=2)
        if share > 8:
            ax.text(left + share / 2, 0, f"{share:.1f} %", ha="center", va="center",
                    color="white", fontsize=11, fontweight="bold")
        left += share
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, 0.75)
    ax.set_yticks([])
    ax.set_xlabel("share of the total temperature rise [%]")
    ax.set_title(f"{total:.1f} K from inlet air to winding mean", fontsize=10.5)
    ax.legend(handles=[Patch(color=c, label=l.replace("\n", " ")) for l, _, c in stages],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, frameon=False,
              fontsize=8.6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    ax.set_axisbelow(True)

    # --- right: the same three numbers, log scale ------------------------
    ax = axes[1]
    labels = ["air bulk\nrise", "air film", "metal"]
    values = [m["dt_mixed_K"], m["superheat_K"], m["metal_drop_K"]]
    colours = [BLUE, ORANGE, COPPER]
    ypos = range(len(values))
    ax.barh(list(ypos), values, color=colours, height=0.5, zorder=2)
    ax.set_xscale("log")
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(1e-2, 3e2)
    for y, v in zip(ypos, values):
        ax.text(v * 1.15, y, f"{v:.3f} K", va="center", fontsize=9)
    ax.set_xlabel("temperature increment [K], log scale")
    ax.set_title("the same increments, logarithmic", fontsize=10.5)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax.text(0.5, -0.42,
            f"wetted surface {m['t_wetted_avg_K'] - 273.15:.2f} °C   "
            f"winding peak {m['t_winding_max_K'] - 273.15:.2f} °C   "
            f"air in {m['t_in_K'] - 273.15:.2f} °C",
            transform=ax.transAxes, ha="center", fontsize=8.6, family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="#f4f6f8",
                      edgecolor="#d5dbe1"))

    fig.suptitle("Where the temperature rise actually happens", fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(FIG / "thermal_path.png")
    plt.close(fig)


def fig_validation(d: dict) -> None:
    """Predicted vs measured on the two quantities the design point fixed."""
    m = d["measured"]
    p = d["prediction"]

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.4), dpi=160)

    ax = axes[0]
    labels = ["pre-run\nprediction", "CFD"]
    vals = [m["housing_over_inlet_predicted_K"], m["housing_over_inlet_K"]]
    bars = ax.bar(labels, vals, color=[GREY, ORANGE], width=0.5, zorder=2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}", ha="center",
                va="bottom", fontsize=10)
    ax.set_ylabel("wetted surface - inlet air [K]")
    ax.set_title(f"Surface superheat   {m['housing_vs_prediction_pct']:+.1f} %",
                 fontsize=10.5)
    ax.set_ylim(0, max(vals) * 1.24)
    style(ax)

    ax = axes[1]
    vals = [p["h_W_m2K"], m["h_effective_W_m2K"]]
    bars = ax.bar(labels, vals, color=[GREY, BLUE], width=0.5, zorder=2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}", ha="center",
                va="bottom", fontsize=10)
    ax.set_ylabel("h [W/m²K]")
    ax.set_title(f"Effective h   {m['h_vs_prediction_pct']:+.1f} %", fontsize=10.5)
    ax.set_ylim(0, max(vals) * 1.24)
    style(ax)

    fig.suptitle("Conjugate heat transfer: prediction fixed before the mesh vs CFD",
                 fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "validation_prediction.png")
    plt.close(fig)


def fig_energy(d: dict) -> None:
    """The area-average / mixed-mean gap, and why it is not an imbalance."""
    m = d["measured"]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), dpi=160)

    ax = axes[0]
    vals = [m["dt_area_K"], m["dt_mixed_K"]]
    bars = ax.bar(["area average\n(as measured)", "mixed mean\n(implied by 80 W)"],
                  vals, color=[ORANGE, GREEN], width=0.5, zorder=2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.4f} K", ha="center",
                va="bottom", fontsize=10)
    ax.set_ylabel("outlet - inlet air temperature [K]")
    ax.set_title(f"Area average overstates the bulk rise by "
                 f"{m['area_vs_mixed_bias_pct']:+.1f} %", fontsize=10.5)
    ax.set_ylim(0, max(vals) * 1.26)
    style(ax)

    ax = axes[1]
    vals = [80.0, m["q_area_W"], 80.0]
    labels = ["imposed\n(power source)", "area average\n(naive)", "mixed mean\n(by definition)"]
    bars = ax.bar(labels, vals, color=[GREY, ORANGE, GREEN], width=0.5, zorder=2)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f} W", ha="center",
                va="bottom", fontsize=9.5)
    ax.set_ylabel("heat carried by the air [W]")
    ax.set_title("Same field, three ways of averaging it", fontsize=10.5)
    ax.set_ylim(0, max(vals) * 1.24)
    style(ax)

    fig.suptitle("The outlet area average is not the bulk temperature", fontsize=11.5,
                 y=0.99)
    fig.text(0.5, 0.005,
             "Not an energy imbalance: the finite-volume solve conserves energy by construction. "
             "The slow warm wake covers more area than it does mass flow.",
             ha="center", fontsize=8.2, color=GREY)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(FIG / "energy_balance.png")
    plt.close(fig)


def fig_convergence(d: dict) -> None:
    p = RESULTS / "convergence.json"
    if not p.exists():
        print("  (no convergence.json - skipping convergence.png)")
        return
    c = json.loads(p.read_text(encoding="utf-8"))

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), dpi=160)

    ax = axes[0]
    for name, series in c["monitors"].items():
        xs = [pt[0] for pt in series]
        final = series[-1][1]
        # Absolute deviation from the converged value, relative to that value,
        # so three monitors in K, Pa and m/s can share one axis.
        ys = []
        for pt in series:
            d = abs(pt[1] - final) / max(abs(final), 1e-12)
            # Values are stored to 6 decimals, so a monitor that has stopped
            # moving reads exactly equal to its final value.  Plot those at the
            # readout floor rather than dropping them, and mark the floor.
            ys.append(max(d, 1e-8))
        ax.plot(xs, ys, lw=1.1, label=name)
    ax.axhline(1e-8, color=GREY, ls=":", lw=1.0, zorder=1)
    ax.text(1000, 1.15e-8, "6-decimal readout floor", ha="right", va="bottom",
            fontsize=7.2, color=GREY)
    ax.set_yscale("log")
    ax.set_xlabel("iteration")
    ax.set_ylabel("|x - x_final| / |x_final|")
    ax.set_title("Surface monitors: relative drift", fontsize=10.5)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, zorder=0)

    ax = axes[1]
    for name, series in c["residuals"].items():
        xs = [pt[0] for pt in series]
        # A residual that has fallen below the readout resolution (1e-6) is
        # stored as 0.0; break the line there rather than invent a floor.
        ys = [float("nan") if pt[1] <= 0 else pt[1] for pt in series]
        ax.plot(xs, ys, lw=1.0, label=name)
    ax.set_yscale("log")
    ax.set_xlabel("iteration")
    ax.set_ylabel("scaled residual")
    ax.set_title("Residuals", fontsize=10.5)
    ax.legend(fontsize=7.5, ncol=2, frameon=False)
    ax.grid(alpha=0.25, zorder=0)
    ax.text(0.98, 0.03, "ω falls below the 1e-6 readout\nresolution and the line ends",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.2, color=GREY)

    fig.suptitle(f"Convergence: {c['iterations']} iterations, "
                 f"steady-state pseudo-time", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "convergence.png")
    plt.close(fig)


def main() -> None:
    d = load()
    FIG.mkdir(exist_ok=True)
    fig_thermal_path(d)
    fig_validation(d)
    fig_energy(d)
    fig_convergence(d)
    print("wrote figures/thermal_path.png, validation_prediction.png, "
          "energy_balance.png, convergence.png")


if __name__ == "__main__":
    main()
