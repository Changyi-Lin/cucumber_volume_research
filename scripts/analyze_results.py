#!/usr/bin/env python3
"""Summarize benchmark CSVs and repeated-run distributions."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import statistics
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def percentile(values: list[float], q: float) -> float:
    return float(np.percentile(np.asarray(values, dtype=float), q))


def distribution(values: list[float]) -> dict[str, float]:
    return {
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "std": statistics.pstdev(values),
        "min": min(values),
        "max": max(values),
        "p95": percentile(values, 95),
    }


def repeat_summary() -> dict[str, dict[str, dict[str, float]]]:
    sources = {
        "balanced": RESULTS / "repeat_balanced.csv",
        "fastest": RESULTS / "repeat_fastest.csv",
        "quality": RESULTS / "repeat_quality.csv",
    }
    metrics = ["downsample_ms", "alpha_wrap_ms", "component_cleanup_ms",
               "volume_ms", "reconstruction_only_ms", "online_core_ms",
               "validation_ms", "total_ms"]
    summary: dict[str, dict[str, dict[str, float]]] = {}
    output_rows = []
    for name, path in sources.items():
        rows = read_rows(path)
        summary[name] = {}
        for metric in metrics:
            stats = distribution([number(row, metric) for row in rows])
            summary[name][metric] = stats
            output_rows.append({"candidate": name, "metric": metric, **stats})
    with (RESULTS / "repeat_statistics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    return summary


def find_row(rows: list[dict[str, str]], voxel: float, alpha: float,
             offset: float, valid_only: bool = True) -> dict[str, str]:
    candidates = [r for r in rows
                  if abs(number(r, "voxel_size_mm") - voxel) < 1e-8
                  and abs(number(r, "alpha_mm") - alpha) < 1e-8
                  and abs(number(r, "offset_mm") - offset) < 1e-8
                  and (not valid_only or r["valid"] == "1")]
    if not candidates:
        raise KeyError((voxel, alpha, offset))
    return candidates[-1]


def main() -> None:
    rows = read_rows(RESULTS / "alpha_wrap_benchmark.csv")
    point_rows = sorted(read_rows(RESULTS / "point_count_runtime.csv"),
                        key=lambda r: number(r, "input_points"))
    repeat = repeat_summary()
    valid = [r for r in rows if r["valid"] == "1"]
    invalid = [r for r in rows if r["valid"] != "1"]

    selected = {
        "fastest": find_row(rows, 3.0, 40.0, 0.75),
        "balanced": find_row(rows, 1.0, 6.25, 0.1),
        "quality": find_row(rows, 0.0, 5.0, 0.1),
    }
    plateau_rows = [find_row(rows, voxel, 6.25, 0.1)
                    for voxel in [1.0, 1.25, 1.5, 1.75, 2.0]]
    plateau_volumes = [number(r, "volume_ml") for r in plateau_rows]
    plateau = {
        "configuration": "alpha=6.25 mm, offset=0.10 mm, voxel=1.0..2.0 mm",
        "volumes_ml": plateau_volumes,
        "mean_ml": statistics.fmean(plateau_volumes),
        "std_ml": statistics.pstdev(plateau_volumes),
        "cv_percent": 100 * statistics.pstdev(plateau_volumes) / statistics.fmean(plateau_volumes),
        "range_ml": max(plateau_volumes) - min(plateau_volumes),
        "relative_range_percent": 100 * (max(plateau_volumes) - min(plateau_volumes)) /
                                  statistics.fmean(plateau_volumes),
    }

    below = [r for r in point_rows if number(r, "online_core_ms") < 500]
    above = [r for r in point_rows if number(r, "online_core_ms") >= 500]
    point_threshold = {
        "largest_measured_below_500_points": int(number(below[-1], "input_points")),
        "largest_measured_below_500_ms": number(below[-1], "online_core_ms"),
        "first_measured_above_500_points": int(number(above[0], "input_points")),
        "first_measured_above_500_ms": number(above[0], "online_core_ms"),
        "recommended_margin_points": 750000,
    }

    compact_selected = {}
    for name, row in selected.items():
        compact_selected[name] = {
            "voxel_size_mm": number(row, "voxel_size_mm"),
            "input_points": int(number(row, "input_points")),
            "alpha_mm": number(row, "alpha_mm"),
            "offset_mm": number(row, "offset_mm"),
            "volume_ml": number(row, "volume_ml"),
            "vertices": int(number(row, "vertices")),
            "faces": int(number(row, "faces")),
            "surface_area_mm2": number(row, "surface_area_mm2"),
            "watertight": bool(int(row["watertight"])),
            "manifold": bool(int(row["manifold"])),
            "two_sided_wrap": bool(int(row["two_sided_wrap"])),
        }

    summary = {
        "experiment_rows": len(rows),
        "valid_rows": len(valid),
        "invalid_rows": len(invalid),
        "two_sided_rows": sum(r["two_sided_wrap"] == "1" for r in rows),
        "selected": compact_selected,
        "repeat_statistics": repeat,
        "volume_plateau": plateau,
        "point_count_threshold": point_threshold,
        "original_alpha15_offset03": {
            key: number(point_rows[-1], key)
            for key in ["alpha_wrap_ms", "volume_ms", "online_core_ms", "volume_ml"]
        },
        "two_sided_cases": [
            {key: r[key] for key in ["voxel_size_mm", "alpha_mm", "offset_mm",
                                      "volume_ml", "failure_reason"]}
            for r in invalid if r["two_sided_wrap"] == "1"
        ],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
