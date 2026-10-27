"""
Render docs/report.md from docs/report_template.md and results/summary.json.

The template is the source of truth for the prose; every number in it is a
placeholder of the form ``{{a.b.c}}`` looked up in the summary JSON.  A
missing key is a hard error, so the report can never silently contain a
stale or blank figure.

Formatting is applied by suffix convention:

    {{measured.h_effective_W_m2K}}        plain repr
    {{measured.h_effective_W_m2K|1}}      one decimal place
    {{measured.h_effective_W_m2K|0}}      integer
    {{measured.h_effective_W_m2K|pct}}    signed percentage, one decimal
    {{measured.h_effective_W_m2K|K}}      kelvin -> degrees Celsius
    {{mesh.cells|num}}                    thousands separators

Usage
-----
    python fill_report.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOCS = ROOT / "docs"
TEMPLATE = DOCS / "report_template.md"
SUMMARY = ROOT / "results" / "summary.json"
OUT = DOCS / "report.md"

PLACEHOLDER = re.compile(r"\{\{([^}]+)\}\}")


def resolve(data: dict, path: str):
    node = data
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(path)
        node = node[part]
    return node


def fmt(value, spec: str | None) -> str:
    if value is None:
        return "n/a"
    if spec is None:
        if isinstance(value, float):
            return f"{value:g}"
        return str(value)
    if spec == "pct":
        return f"{value:+.1f}"
    if spec == "K":
        return f"{value - 273.15:.1f}"
    if spec == "K2":
        return f"{value - 273.15:.2f}"
    if spec.isdigit():
        return f"{value:.{spec}f}"
    if spec == "num":
        return f"{value:,.0f}"
    if spec == "e":
        return f"{value:.2e}"
    raise ValueError(f"unknown format spec {spec!r}")


def main() -> None:
    if not TEMPLATE.exists():
        raise SystemExit(f"missing {TEMPLATE}")
    if not SUMMARY.exists():
        raise SystemExit(f"missing {SUMMARY} - run derive.py first")

    data = json.loads(SUMMARY.read_text(encoding="utf-8"))
    text = TEMPLATE.read_text(encoding="utf-8")

    missing: list[str] = []

    def repl(m: re.Match) -> str:
        expr = m.group(1)
        path, _, spec = expr.partition("|")
        try:
            return fmt(resolve(data, path.strip()), spec.strip() or None)
        except KeyError:
            missing.append(path.strip())
            return m.group(0)

    rendered = PLACEHOLDER.sub(repl, text)
    if missing:
        uniq = sorted(set(missing))
        raise SystemExit("template refers to missing keys:\n  " + "\n  ".join(uniq))

    OUT.write_text(rendered, encoding="utf-8", newline="\n")
    print(f"wrote {OUT}  ({len(rendered):,} chars)")


if __name__ == "__main__":
    main()
