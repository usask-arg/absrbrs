"""``python -m faddeeva_prototype report``: measure the baselines, search, and write the report."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from faddeeva_prototype import occupancy, reference, search

DOCS = reference.REPO / "docs" / "validation"
DEFAULT_JSON = DOCS / "data" / "faddeeva_prototype.json"
DEFAULT_TABLE = DOCS / "faddeeva-prototype-table.md"


def _fmt(v: float) -> str:
    return f"{v:.1e}"


def table_markdown(
    baselines: dict[str, search.Result],
    best: search.Result | None,
    holdout: dict[str, tuple[float, float]],
) -> str:
    lines = [
        "| Algorithm | Max value error | Max derivative error | Holdout value error "
        "| Holdout derivative error | Expected cost |",
        "|---|---|---|---|---|---|",
    ]
    rows = [*baselines.items(), *([("**Selected composite**", best)] if best else [])]
    for name, r in rows:
        hv, hd = holdout[name]
        lines.append(
            f"| {name} | {_fmt(r.max_value_error)} | {_fmt(r.max_derivative_error)} "
            f"| {_fmt(hv)} | {_fmt(hd)} | {r.expected_cost:.0f} |"
        )
    if best:
        lines += [
            "",
            "Regions of the selected composite (first match wins):",
            "",
            "| Region | Method | Max value error | Max derivative error | Workload share |",
            "|---|---|---|---|---|",
        ]
        for region, errors, share in zip(
            best.regions, best.region_errors, best.occupancy, strict=True
        ):
            method = f"{region.method[0]} {region.method[1]}"
            name = region.name.replace("|", r"\|")  # names such as "CF |z| in [8, 10)"
            lines.append(
                f"| {name} | {method} | {_fmt(errors[0])} | {_fmt(errors[1])} | {share:.2%} |"
            )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="faddeeva_prototype")
    sub = parser.add_subparsers(dest="command", required=True)
    rep = sub.add_parser("report", help="measure, search and write the report")
    rep.add_argument("--json", type=Path, default=DEFAULT_JSON)
    rep.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    args = parser.parse_args(argv)

    zc, wc, dwc = reference.load_committed()
    zd, wd, dwd = reference.dense()
    zh, wh, dwh = reference.holdout()
    z, w, dw = (np.concatenate(p) for p in ((zc, zd), (wc, wd), (dwc, dwd)))
    table = search.ErrorTable(z, w, dw)
    workload = occupancy.workload_points()

    baselines = {
        name: search.measure(regions, table, workload) for name, regions in search.BASELINES.items()
    }
    holdout = {
        name: search.holdout_errors(search.BASELINES[name], zh, wh, dwh) for name in baselines
    }
    results = search.search(table, workload)
    # The selected composite is the cheapest that also meets the targets, with the margin, on
    # the random holdout: the grid alone could miss an error peak between its points.
    rejected = 0
    best = None
    for candidate in results:
        errors = search.holdout_errors(candidate.regions, zh, wh, dwh)
        if search.meets(errors):
            best = candidate
            holdout["**Selected composite**"] = errors
            break
        rejected += 1
    report = {
        "targets": {
            "value": search.TARGET_VALUE,
            "derivative": search.TARGET_DERIVATIVE,
            "margin": search.MARGIN,
        },
        "evaluation_points": {
            "committed": int(zc.size),
            "dense": int(zd.size),
            "holdout": int(zh.size),
        },
        "workload_points": int(workload.size),
        "baselines": {
            name: {**r.to_json(), "holdout": [search.json_float(e) for e in holdout[name]]}
            for name, r in baselines.items()
        },
        "valid_composites": len(results),
        "rejected_on_holdout": rejected,
        "selected": (
            {
                **best.to_json(),
                "holdout": [search.json_float(e) for e in holdout["**Selected composite**"]],
            }
            if best
            else None
        ),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    markdown = table_markdown(baselines, best, holdout)
    args.table.write_text(markdown, encoding="utf-8")
    print(markdown, file=sys.stderr)
    if best is None:
        print("No composite meets the targets: see the baselines in the report.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
