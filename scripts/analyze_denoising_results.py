#!/usr/bin/env python3
"""Create deterministic summary tables from the denoising benchmark CSVs."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
METHODS = ("NONE", "SOR", "ROR", "SOR_ROR")
VOXELS = (1.0, 1.5, 2.0, 2.5)
ALPHAS = (6.25, 7.5, 10.0)
OFFSETS = (0.1, 0.2, 0.3)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in ("", None) else float("nan")


def close(a: float, b: float) -> bool:
    return abs(a - b) < 1e-8


def selected_config(row: dict[str, str]) -> bool:
    method = row["denoise_method"]
    voxel = number(row, "voxel_size_mm")
    if method == "NONE":
        return close(number(row, "end_protection"), 0.0)
    if method == "SOR":
        return (int(number(row, "sor_k")) == 30 and
                close(number(row, "sor_std_ratio"), 3.0) and
                close(number(row, "end_protection"), 0.05))
    if method == "ROR":
        return (close(number(row, "ror_radius_mm"), 2.5 * voxel) and
                int(number(row, "ror_min_neighbors")) == 3 and
                close(number(row, "end_protection"), 0.05))
    if method == "SOR_ROR":
        return (int(number(row, "sor_k")) == 30 and
                close(number(row, "sor_std_ratio"), 3.0) and
                close(number(row, "ror_radius_mm"), 2.5 * voxel) and
                int(number(row, "ror_min_neighbors")) == 3 and
                close(number(row, "end_protection"), 0.05))
    return False


def clean_selected(row: dict[str, str]) -> bool:
    return (selected_config(row) and
            close(number(row, "injected_outlier_percent"), 0.0) and
            int(number(row, "clustering")) == 0)


def deduplicate(rows: list[dict[str, str]], keys: tuple[str, ...]) -> list[dict[str, str]]:
    latest: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        latest[key] = row
    return list(latest.values())


def stats(values: list[float]) -> dict[str, float]:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return {name: float("nan") for name in
                ("mean", "std", "min", "max", "range", "cv_percent", "p95")}
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=1)) if data.size > 1 else 0.0
    return {
        "mean": mean,
        "std": std,
        "min": float(np.min(data)),
        "max": float(np.max(data)),
        "range": float(np.ptp(data)),
        "cv_percent": 100.0 * std / mean if mean else float("nan"),
        "p95": float(np.percentile(data, 95)),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def comparison_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates = [row for row in rows if clean_selected(row)
                  and number(row, "voxel_size_mm") in VOXELS
                  and number(row, "alpha_mm") in ALPHAS
                  and number(row, "offset_mm") in OFFSETS]
    result = deduplicate(candidates,
                         ("voxel_size_mm", "denoise_method", "alpha_mm", "offset_mm"))
    expected = len(VOXELS) * len(METHODS) * len(ALPHAS) * len(OFFSETS)
    if len(result) != expected:
        raise RuntimeError(f"Expected {expected} comparison rows, found {len(result)}")
    return result


def build_method_summary(rows: list[dict[str, str]], repeats: list[dict[str, str]]) -> None:
    comparison = comparison_rows(rows)
    output: list[dict[str, object]] = []
    for method in METHODS:
        grid = [row for row in comparison if row["denoise_method"] == method]
        representative = [row for row in grid if close(number(row, "voxel_size_mm"), 1.5)]
        repeated = [row for row in repeats if row["denoise_method"] == method]
        volume = stats([number(row, "volume_ml") for row in grid])
        representative_volume = stats([number(row, "volume_ml") for row in representative])
        within_voxel_cv = []
        for voxel in VOXELS:
            values = [number(row, "volume_ml") for row in grid
                      if close(number(row, "voxel_size_mm"), voxel)]
            within_voxel_cv.append(stats(values)["cv_percent"])
        repeat_denoise = stats([number(row, "denoise_total_ms") for row in repeated])
        repeat_wrap = stats([number(row, "alpha_wrap_ms") for row in repeated])
        repeat_deploy = stats([number(row, "deployment_total_ms") for row in repeated])
        repeat_volume = stats([number(row, "volume_ml") for row in repeated])
        sample = repeated[-1]
        output.append({
            "method": method,
            "grid_cases": len(grid),
            "grid_valid_percent": 100.0 * sum(int(number(r, "valid")) for r in grid) / len(grid),
            "voxel_size_mm_representative": 1.5,
            "remaining_points": int(number(sample, "points_after")),
            "removed_percent": number(sample, "removed_percent"),
            "top_removed_percent": number(sample, "top_removed_percent"),
            "middle_removed_percent": number(sample, "middle_removed_percent"),
            "bottom_removed_percent": number(sample, "bottom_removed_percent"),
            "length_change_mm": number(sample, "length_change_mm"),
            "repeat_n": len(repeated),
            "denoise_ms_mean": repeat_denoise["mean"],
            "denoise_ms_std": repeat_denoise["std"],
            "alpha_wrap_ms_mean": repeat_wrap["mean"],
            "alpha_wrap_ms_std": repeat_wrap["std"],
            "deployment_total_ms_mean": repeat_deploy["mean"],
            "deployment_total_ms_std": repeat_deploy["std"],
            "deployment_total_ms_p95": repeat_deploy["p95"],
            "volume_ml_fixed_mean": repeat_volume["mean"],
            "volume_ml_fixed_std": repeat_volume["std"],
            "volume_ml_grid_mean": volume["mean"],
            "volume_ml_grid_std": volume["std"],
            "volume_ml_grid_range": volume["range"],
            "volume_cv_grid_percent": volume["cv_percent"],
            "volume_cv_mean_within_voxel_percent": float(np.mean(within_voxel_cv)),
            "volume_cv_v1p5_neighborhood_percent": representative_volume["cv_percent"],
            "vertices_fixed": int(number(sample, "vertices")),
            "faces_fixed": int(number(sample, "faces")),
            "watertight": int(number(sample, "watertight")),
            "manifold": int(number(sample, "manifold")),
            "two_sided_wrap": int(number(sample, "two_sided_wrap")),
        })
    write_csv(RESULTS / "denoising_method_summary.csv", output)


def build_repeat_summary(repeats: list[dict[str, str]]) -> None:
    output: list[dict[str, object]] = []
    metrics = ("voxel_ms", "denoise_total_ms", "alpha_wrap_ms", "volume_ms",
               "deployment_total_ms", "research_total_ms", "volume_ml")
    for method in METHODS:
        subset = [row for row in repeats if row["denoise_method"] == method]
        record: dict[str, object] = {"method": method, "runs": len(subset)}
        for metric in metrics:
            summary = stats([number(row, metric) for row in subset])
            record[f"{metric}_mean"] = summary["mean"]
            record[f"{metric}_std"] = summary["std"]
            record[f"{metric}_p95"] = summary["p95"]
        output.append(record)
    write_csv(RESULTS / "denoising_repeatability_summary.csv", output)


def build_stress(rows: list[dict[str, str]]) -> None:
    baselines = {}
    for row in rows:
        if (clean_selected(row) and close(number(row, "voxel_size_mm"), 1.5)
                and close(number(row, "alpha_mm"), 6.25)
                and close(number(row, "offset_mm"), 0.1)):
            baselines[row["denoise_method"]] = row
    stress_rows = [row for row in rows if selected_config(row)
                   and number(row, "injected_outlier_percent") > 0
                   and close(number(row, "voxel_size_mm"), 1.5)]
    if len(stress_rows) != 60:
        raise RuntimeError(f"Expected 60 contaminated stress rows, found {len(stress_rows)}")
    detailed: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in stress_rows:
        method = row["denoise_method"]
        baseline_volume = number(baselines[method], "volume_ml")
        volume = number(row, "volume_ml")
        record = {
            "method": method,
            "outlier_percent": number(row, "injected_outlier_percent"),
            "outlier_distance_mm": number(row, "injected_outlier_distance_mm"),
            "points_before": int(number(row, "points_before")),
            "removed_points": int(number(row, "removed_points")),
            "removed_percent": number(row, "removed_percent"),
            "volume_ml": volume,
            "baseline_volume_ml": baseline_volume,
            "volume_delta_ml": volume - baseline_volume,
            "volume_delta_percent": 100.0 * (volume - baseline_volume) / baseline_volume,
            "alpha_wrap_ms": number(row, "alpha_wrap_ms"),
            "deployment_total_ms": number(row, "deployment_total_ms"),
            "watertight": int(number(row, "watertight")),
            "two_sided_wrap": int(number(row, "two_sided_wrap")),
            "valid": int(number(row, "valid")),
            "failure_reason": row["failure_reason"],
        }
        detailed.append(record)
        grouped[method].append(record)
    write_csv(RESULTS / "outlier_stress_analysis.csv", detailed)

    summary: list[dict[str, object]] = []
    for method in METHODS:
        subset = grouped[method]
        harsh = next(r for r in subset
                     if close(float(r["outlier_percent"]), 1.0)
                     and close(float(r["outlier_distance_mm"]), 20.0))
        abs_delta = [abs(float(r["volume_delta_percent"])) for r in subset]
        summary.append({
            "method": method,
            "stress_cases": len(subset),
            "valid_percent": 100.0 * sum(int(r["valid"]) for r in subset) / len(subset),
            "mean_abs_volume_delta_percent": float(np.mean(abs_delta)),
            "worst_abs_volume_delta_percent": float(np.max(abs_delta)),
            "harsh_1pct_20mm_removed_points": harsh["removed_points"],
            "harsh_1pct_20mm_volume_ml": harsh["volume_ml"],
            "harsh_1pct_20mm_volume_delta_percent": harsh["volume_delta_percent"],
            "harsh_1pct_20mm_deployment_total_ms": harsh["deployment_total_ms"],
            "harsh_1pct_20mm_valid": harsh["valid"],
            "harsh_1pct_20mm_failure_reason": harsh["failure_reason"],
        })
    write_csv(RESULTS / "outlier_stress_summary.csv", summary)


def build_end_protection(rows: list[dict[str, str]]) -> None:
    latest: dict[tuple[str, float], dict[str, str]] = {}
    for row in rows:
        method = row["denoise_method"]
        end = number(row, "end_protection")
        if (method not in ("SOR", "ROR", "SOR_ROR") or
                end not in (0.0, 0.05, 0.075, 0.1) or
                not close(number(row, "voxel_size_mm"), 1.5) or
                not close(number(row, "alpha_mm"), 6.25) or
                not close(number(row, "offset_mm"), 0.1) or
                number(row, "injected_outlier_percent") != 0 or
                int(number(row, "clustering")) != 0):
            continue
        if method in ("SOR", "SOR_ROR") and not (
                int(number(row, "sor_k")) == 30 and close(number(row, "sor_std_ratio"), 3.0)):
            continue
        if method in ("ROR", "SOR_ROR") and not (
                close(number(row, "ror_radius_mm"), 3.75) and
                int(number(row, "ror_min_neighbors")) == 3):
            continue
        latest[(method, end)] = row
    output = []
    for (method, end), row in sorted(latest.items()):
        output.append({key: row[key] for key in
                       ("denoise_method", "end_protection", "points_before", "points_after",
                        "removed_percent", "top_removed_percent", "middle_removed_percent",
                        "bottom_removed_percent", "length_change_mm", "length_change_percent",
                        "possible_over_denoise", "denoise_total_ms", "alpha_wrap_ms",
                        "deployment_total_ms", "volume_ml", "valid", "failure_reason")})
    if len(output) != 12:
        raise RuntimeError(f"Expected 12 end-protection rows, found {len(output)}")
    write_csv(RESULTS / "end_protection_summary.csv", output)


def build_clustering(rows: list[dict[str, str]]) -> None:
    output = []
    for row in rows:
        if int(number(row, "clustering")) != 1:
            continue
        output.append({key: row[key] for key in
                       ("end_protection", "cluster_radius_mm", "number_of_clusters",
                        "largest_cluster_points", "second_cluster_points", "removed_cluster_points",
                        "removed_percent", "length_change_mm", "clustering_ms", "volume_ml",
                        "valid", "failure_reason")})
    write_csv(RESULTS / "clustering_summary.csv", output)


def main() -> None:
    rows = read_rows(RESULTS / "denoising_benchmark.csv")
    repeats = read_rows(RESULTS / "denoising_repeatability.csv")
    build_method_summary(rows, repeats)
    build_repeat_summary(repeats)
    build_stress(rows)
    build_end_protection(rows)
    build_clustering(rows)
    print(f"Analyzed {len(rows)} benchmark rows and {len(repeats)} repeat rows.")


if __name__ == "__main__":
    main()
