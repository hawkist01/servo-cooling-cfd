"""
Geometry preview for the servo motor CHT case.

Two panels, both drawn to scale from the exported meshes:

  left   the flow domain on the z = 0 plane: air, the motor inside it, and the
         inlet / outlet boundaries
  right  the motor itself on the same plane, dimensioned

The housing is closed at both ends and the winding is fully enclosed, so there
is no wetted winding surface: the copper never meets the air. The section panel
shows exactly one conduction interface and one convection surface, and nothing
else.

Run:  python preview.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
OUT = HERE / "geometry_preview.png"

COPPER = "#c1440e"
ALUM = "#8aa6c0"
AIRFILL = "#dbe8f4"
BLUE = "#1f6fb2"
EDGE = "#31404d"

# dimensions, mm -- must match make_geometry.py
R_WIND = 26.0
L_WIND = 112.0
R_SLEEVE_OUT = 32.0
R_FIN_TIP = 42.0
T_END = 2.0
L_HOUSING = L_WIND + 2.0 * T_END
X_UP, X_DOWN = 150.0, 450.0
HALF_W = 150.0


def section_polygons(mesh: trimesh.Trimesh):
    """Closed 2D polygons of the z = 0 section, in mm.

    to_2D() returns coordinates in a plane basis of its own choosing, which for
    this mesh came back reflected in x and shifted by 257 mm -- the section
    looked plausible but was not in model coordinates. Each polygon is mapped
    back through the returned transform and then read off as (x, y), so the
    numbers are the model's own.
    """
    sec = mesh.section(plane_origin=[0, 0, 0], plane_normal=[0, 0, 1])
    planar, to_3d = sec.to_2D()

    def to_xy(xy: np.ndarray) -> np.ndarray:
        pts = np.column_stack([xy, np.zeros(len(xy))])
        return trimesh.transform_points(pts, to_3d)[:, :2] * 1e3

    out = []
    for poly in planar.polygons_full:
        out.append((
            to_xy(np.asarray(poly.exterior.xy).T),
            [to_xy(np.asarray(h.xy).T) for h in poly.interiors],
        ))
    return out


def fill_section(ax, mesh: trimesh.Trimesh, face: str, hole_face: str | None = None,
                 zorder: int = 2, lw: float = 0.9) -> None:
    for exterior, holes in section_polygons(mesh):
        ax.fill(exterior[:, 0], exterior[:, 1], facecolor=face,
                edgecolor=EDGE, lw=lw, zorder=zorder)
        # Holes are separate patches drawn ON TOP of the exterior, not behind
        # it: filling the exterior fills the hole too, so a hole drawn first is
        # simply covered again.
        for hole in holes:
            ax.fill(hole[:, 0], hole[:, 1],
                    facecolor=hole_face if hole_face else "white",
                    edgecolor=EDGE, lw=lw, zorder=zorder + 0.5)


def main() -> None:
    winding = trimesh.load(HERE / "motor_winding.stl")
    housing = trimesh.load(HERE / "motor_housing.stl")
    air = trimesh.load(HERE / "air_domain.stl")

    fig = plt.figure(figsize=(14.0, 5.8))

    # ---- flow domain ----------------------------------------------------
    ax1 = fig.add_subplot(1, 2, 1)

    # The air section is the box with the motor cut out, so its holes are the
    # motor's own silhouette -- filling them with the housing colour draws the
    # motor without redrawing it by hand.
    fill_section(ax1, air, AIRFILL, hole_face=ALUM, zorder=1)
    fill_section(ax1, winding, COPPER, zorder=3)

    for y in (-120, -70, 0, 70, 120):
        ax1.annotate("", xy=(-122, y), xytext=(-150, y),
                     arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.2),
                     zorder=4)
    ax1.plot([-X_UP, -X_UP], [-HALF_W, HALF_W], color=BLUE, lw=1.6, zorder=4)
    ax1.text(-X_UP, HALF_W + 8, "inlet\n3 m/s, 25 C", ha="left", va="bottom",
             fontsize=8.5, color=BLUE)
    ax1.plot([X_DOWN, X_DOWN], [-HALF_W, HALF_W], color=BLUE, lw=1.6, zorder=4)
    ax1.text(X_DOWN, HALF_W + 8, "outlet\n0 Pa", ha="right", va="bottom",
             fontsize=8.5, color=BLUE)

    ax1.annotate("", xy=(-X_UP, -162), xytext=(X_DOWN, -162),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax1.text(0.5 * (X_DOWN - X_UP), -170, "600 mm", ha="center", va="top",
             fontsize=8, color="0.25")
    ax1.annotate("", xy=(-172, -HALF_W), xytext=(-172, HALF_W),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax1.text(-178, 0, "300 mm", ha="right", va="center", fontsize=8,
             color="0.25", rotation=90)

    ax1.set_xlim(-200, 480)
    ax1.set_ylim(-192, 196)
    ax1.set_axis_off()
    ax1.set_title("flow domain, mid-plane (z = 0)", fontsize=10.5, pad=2)

    # ---- motor section ---------------------------------------------------
    ax2 = fig.add_subplot(1, 2, 2)

    fill_section(ax2, housing, ALUM, zorder=2)
    fill_section(ax2, winding, COPPER, zorder=3)

    # fin roots: the fin and the sleeve wall are the same material and the same
    # colour, so without this hairline the 10 mm fin is indistinguishable from
    # the 6 mm wall it sits on
    for sgn in (1, -1):
        ax2.plot([-L_WIND / 2, L_WIND / 2], [sgn * R_SLEEVE_OUT] * 2,
                 color=EDGE, lw=0.55, zorder=4)

    # dimension annotations, in mm (this panel is scaled to mm above)
    ax2.annotate("", xy=(-L_WIND / 2, -47), xytext=(L_WIND / 2, -47),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax2.text(0, -49.5, f"{L_WIND:.0f} mm winding", ha="center", va="top",
             fontsize=8, color="0.25")

    ax2.annotate("", xy=(-L_HOUSING / 2, 46), xytext=(L_HOUSING / 2, 46),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax2.text(0, 47.5, f"{L_HOUSING:.0f} mm housing, closed both ends",
             ha="center", va="bottom", fontsize=8, color="0.25")

    ax2.annotate("", xy=(62, R_SLEEVE_OUT), xytext=(62, R_FIN_TIP),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax2.text(65, 0.5 * (R_SLEEVE_OUT + R_FIN_TIP), "10 mm fin",
             ha="left", va="center", fontsize=8, color="0.25")

    ax2.annotate("", xy=(-66, 0), xytext=(-66, R_WIND),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax2.text(-69, R_WIND / 2, "r = 26", ha="right", va="center",
             fontsize=8, color="0.25")

    # heat source: labelled on the region itself, in white, so it cannot be
    # confused with a callout for one of the surfaces
    ax2.text(0, 4, "80 W generated\nin the winding", color="white", ha="center",
             va="center", fontsize=9.5, fontweight="bold", zorder=5)

    # conduction interface: arrow into the bore wall
    ax2.annotate("", xy=(22, R_WIND), xytext=(62, 10),
                 arrowprops=dict(arrowstyle="->", color=EDGE, lw=0.9))
    ax2.text(64, 8, "heat conducts\nthrough this\ninterface", ha="left",
             va="center", fontsize=8, color=EDGE)

    ax2.annotate("12 axial fins", xy=(38, 38), xytext=(74, 56),
                 ha="center", fontsize=8.5, color="0.2",
                 arrowprops=dict(arrowstyle="->", color="0.4", lw=0.9))

    ax2.set_xlim(-112, 118)
    ax2.set_ylim(-60, 74)
    ax2.set_axis_off()
    ax2.set_title("motor mid-plane section, dimensions in mm",
                  fontsize=10.5, pad=2)

    fig.text(0.5, 0.035,
             "3 m/s air along +x at 25 C.  The winding is fully enclosed by the "
             "housing, so it has no wetted area: every watt crosses the bore "
             "interface and leaves from the outer surface and fins.",
             ha="center", fontsize=8.5, color="0.35")
    fig.legend(handles=[Patch(facecolor=COPPER, edgecolor=EDGE,
                              label="copper winding (80 W source)"),
                        Patch(facecolor=ALUM, edgecolor=EDGE,
                              label="aluminium housing + fins"),
                        Patch(facecolor=AIRFILL, edgecolor=EDGE, label="air")],
               loc="lower center", ncol=3, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.055))

    fig.tight_layout(rect=(0, 0.09, 1, 1))
    # Aspect is set AFTER tight_layout: laying out the figure rewrites the axes
    # box, and an equal aspect applied before it gets stretched again, which
    # silently turns a scale drawing into a misleading one.
    ax1.set_aspect("equal", adjustable="box")
    ax2.set_aspect("equal", adjustable="box")
    fig.canvas.draw()

    for name, ax in [("domain", ax1), ("section", ax2)]:
        b = ax.get_window_extent()
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        px, py = b.width / (x1 - x0), b.height / (y1 - y0)
        if abs(px - py) > 0.01 * px:
            raise SystemExit(
                f"{name} panel is not to scale: {px:.3f} px/mm in x vs "
                f"{py:.3f} px/mm in y"
            )

    fig.savefig(OUT, dpi=170, bbox_inches="tight")
    print(f"Wrote {OUT}  (both panels to scale)")


if __name__ == "__main__":
    main()
