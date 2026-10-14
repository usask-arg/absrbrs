"""``python -m absrbrs_validation faddeeva``: measure both tiers, write the report and table."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from absrbrs_validation import arb, faddeeva

DOCS = Path(__file__).resolve().parents[2] / "docs" / "validation"


def table(results: list[dict]) -> str:
    lines = [
        "| Tier | Points | Max value error | Target | Max derivative error | Target |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['tier']} | {r['points']} | {r['max_value_error']:.1e} | {r['value_target']:.0e} "
            f"| {r['max_derivative_error']:.1e} | {r['derivative_target']:.0e} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="absrbrs_validation")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("faddeeva", help="both Faddeeva tiers against Arb")
    parser.parse_args(argv)
    z = faddeeva.points()
    w, dw = arb.compute(z)
    results = [faddeeva.measure(z, w, dw, tier) for tier in ("fast", "accurate")]
    (DOCS / "data").mkdir(parents=True, exist_ok=True)
    (DOCS / "data" / "faddeeva-accuracy.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    (DOCS / "faddeeva-accuracy-table.md").write_text(table(results), encoding="utf-8")
    print(table(results), file=sys.stderr)
    ok = all(
        r["max_value_error"] <= r["value_target"]
        and r["max_derivative_error"] <= r["derivative_target"]
        for r in results
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
