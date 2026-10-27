"""
Rebuild every derived artefact from results/measurements.json in one command.

    python build.py

This does not run CFD.  It regenerates:
    results/summary.json      (derive.py)
    figures/*.png             (make_figures.py)
    docs/report.md            (fill_report.py, from report_template.md)
    docs/dashboard.html       (make_dashboard.py)

The order matters: derive writes the JSON the other three read.

Any step failing stops the build, so a stale report can never be silently
published alongside fresh numbers.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"

STEPS = [
    ("derive.py", "measurements.json -> energy balance / temperatures / comparison"),
    ("make_figures.py", "figures"),
    ("fill_report.py", "docs/report.md from report_template.md"),
    ("make_dashboard.py", "docs/dashboard.html"),
]


def main() -> int:
    for script, what in STEPS:
        path = ANALYSIS / script
        if not path.exists():
            print(f"missing: {path}", file=sys.stderr)
            return 1
        print(f"== {script:18s} {what}")
        proc = subprocess.run([sys.executable, str(path)], cwd=ANALYSIS)
        if proc.returncode != 0:
            print(f"\n{script} failed (exit {proc.returncode}) - build stopped.",
                  file=sys.stderr)
            return proc.returncode
    print("\nbuild complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
