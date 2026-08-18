#!/usr/bin/env python3
"""Summarize the advanced cleaning, stress, Alpha Wrap, and repeat experiments."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_advanced_cross_sections import read_xyz_ply


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FULL = RESULTS / "full_cleaning_alpha_wrap_benchmark.csv"
REPEAT = RESULTS / "full_cleaning_repeat_best.csv"
FIGURES = RESULTS / "figures" / "advanced_cleaning"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in ("", None) else float("nan")


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def last_match(rows: list[dict[str, str]], predicate) -> dict[str, str]:
    matches = [row for row in rows if predicate(row)]
    if not matches:
        raise RuntimeError("Expected result row was not found")
    return matches[-1]


def main() -> None:
    rows = read_csv(FULL)
    FIGURES.mkdir(parents=True, exist_ok=True)

    geometry_rows: list[dict[str, object]] = []
    for dataset, path in (("RAW", ROOT.parent / "pointcloud_3100.ply"),
                          ("SOR_ROR", RESULTS / "denoised_points" /
                           "sor_ror_remaining.ply")):
        points = read_xyz_ply(path)
        center = np.mean(points, axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(np.cov(points - center, rowvar=False))
        projected = (points - center) @ eigenvectors
        pca_extents = np.max(projected, axis=0) - np.min(projected, axis=0)
        minimum, maximum = np.min(points, axis=0), np.max(points, axis=0)
        geometry_rows.append({
            "dataset": dataset, "point_count": len(points),
            "aabb_min_x": minimum[0], "aabb_min_y": minimum[1], "aabb_min_z": minimum[2],
            "aabb_max_x": maximum[0], "aabb_max_y": maximum[1], "aabb_max_z": maximum[2],
            "aabb_size_x_mm": maximum[0] - minimum[0],
            "aabb_size_y_mm": maximum[1] - minimum[1],
            "aabb_size_z_mm": maximum[2] - minimum[2],
            "pca_extent_small_mm": pca_extents[0],
            "pca_extent_middle_mm": pca_extents[1],
            "pca_main_axis_length_mm": pca_extents[2],
        })
        del points, projected
    write_rows(RESULTS / "advanced_pointcloud_geometry_summary.csv", geometry_rows)

    density_source = read_csv(RESULTS / "advanced_density_analysis.csv")
    density_rows: list[dict[str, object]] = []
    fractions = (5.0, 7.5, 10.0)
    for index, row in enumerate(density_source):
        block = index // 12
        fraction_index = (index % 12) // 4
        density_rows.append({
            "dataset": "RAW" if block == 0 else "SOR_ROR",
            "end_region_percent": fractions[fraction_index],
            "region": row["region"], "points": row["points"],
            "nn_mean_mm": row["mean"], "nn_median_mm": row["median"],
            "nn_std_mm": row["std"], "nn_p5_mm": row["p5"],
            "nn_p25_mm": row["p25"], "nn_p75_mm": row["p75"],
            "nn_p95_mm": row["p95"],
            "density_runtime_ms": row["density_runtime_ms"],
        })
    write_rows(RESULTS / "advanced_density_regional_summary.csv", density_rows)

    residual_pipelines = ("VOXEL", "VOXEL+COMPONENT",
                          "VOXEL+COMPONENT+ADAPTIVE_LOCAL")
    residual: list[dict[str, object]] = []
    for pipeline in residual_pipelines:
        row = last_match(rows, lambda item, p=pipeline:
            item["denoise_pipeline"] == p and number(item, "voxel_size") == 0 and
            number(item, "points_before") == 59781 and
            number(item, "injected_outlier_percent") == 0 and
            number(item, "alpha_mm") == 6.25 and number(item, "offset_mm") == 0.1)
        residual.append({
            "stage": {"VOXEL": "SOR_ROR_BEFORE",
                      "VOXEL+COMPONENT": "AFTER_COMPONENT",
                      "VOXEL+COMPONENT+ADAPTIVE_LOCAL": "AFTER_ADAPTIVE_LOCAL"}[pipeline],
            "remaining_points": int(number(row, "points_after")),
            "total_removed": int(number(row, "removed_points")),
            "component_removed": int(number(row, "removed_component_points")),
            "local_removed": int(number(row, "removed_local_points")),
            "top_removed_percent": number(row, "top_removed_percent"),
            "middle_removed_percent": number(row, "middle_removed_percent"),
            "bottom_removed_percent": number(row, "bottom_removed_percent"),
            "raw_length_change_mm": number(row, "length_change_mm"),
            "confirmed_noise_length_mm": number(row, "confirmed_component_length_change_mm"),
            "unexplained_length_change_mm": number(row, "unexplained_length_change_mm"),
            "cleaning_ms": number(row, "denoise_ms"),
            "alpha_wrap_ms": number(row, "alpha_wrap_ms"),
            "deployment_total_ms": number(row, "deployment_total_ms"),
            "volume_ml": number(row, "volume_ml"),
            "visible_disconnected_noise": "YES" if pipeline == "VOXEL" else "NO",
            "visual_geometry_damage": "NO",
            "valid_after_manual_component_confirmation": int(number(row, "valid")),
        })
    write_rows(RESULTS / "residual_noise_stage_summary.csv", residual)

    requested = (
        "VOXEL", "VOXEL+SOR", "VOXEL+ROR", "VOXEL+SOR+ROR",
        "VOXEL+SOR+ROR+COMPONENT",
        "VOXEL+SOR+ROR+COMPONENT+ADAPTIVE_LOCAL",
        "VOXEL+COMPONENT+ADAPTIVE_LOCAL",
    )
    ablation: list[dict[str, object]] = []
    for voxel in (1.5, 2.0):
        for pipeline in requested:
            candidates = [row for row in rows if
                row["denoise_pipeline"] == pipeline and
                number(row, "voxel_size") == voxel and
                number(row, "alpha_mm") == 6.25 and number(row, "offset_mm") == 0.1 and
                row["injected_outlier_mode"] == "dispersed" and
                number(row, "injected_outlier_percent") == 0]
            if not candidates:
                raise RuntimeError(f"Missing ablation: {voxel}, {pipeline}")
            row = candidates[0]
            ablation.append({
                "voxel_size_mm": voxel,
                "pipeline": pipeline,
                "remaining_points": int(number(row, "points_after")),
                "removed_percent": number(row, "removed_percent"),
                "component_removed": int(number(row, "removed_component_points")),
                "local_removed": int(number(row, "removed_local_points")),
                "end_damage_mm_unexplained": number(row, "unexplained_length_change_mm"),
                "cleaning_ms": number(row, "denoise_ms"),
                "alpha_wrap_ms": number(row, "alpha_wrap_ms"),
                "deployment_total_ms": number(row, "deployment_total_ms"),
                "volume_ml": number(row, "volume_ml"),
                "vertices": int(number(row, "vertices")),
                "faces": int(number(row, "faces")),
                "watertight": int(number(row, "watertight")),
                "manifold": int(number(row, "manifold")),
                "two_sided_wrap": int(number(row, "two_sided_wrap")),
                "valid": int(number(row, "valid")),
            })
    write_rows(RESULTS / "advanced_ablation_summary.csv", ablation)

    # Deduplicate the sweep because final mesh export repeats three candidates.
    alpha_map: dict[tuple[float, float], dict[str, str]] = {}
    for row in rows:
        if (row["denoise_pipeline"] == "VOXEL+COMPONENT" and
                number(row, "voxel_size") == 1.5 and
                number(row, "confirmed_disconnected_noise") == 1 and
                number(row, "injected_outlier_percent") == 0 and
                row["injected_outlier_mode"] == "dispersed"):
            key = (number(row, "alpha_mm"), number(row, "offset_mm"))
            alpha_map.setdefault(key, row)
    alpha_rows: list[dict[str, object]] = []
    for (alpha, offset), row in sorted(alpha_map.items()):
        alpha_rows.append({
            "alpha_mm": alpha, "offset_mm": offset,
            "alpha_wrap_ms": number(row, "alpha_wrap_ms"),
            "deployment_total_ms": number(row, "deployment_total_ms"),
            "volume_ml": number(row, "volume_ml"),
            "vertices": int(number(row, "vertices")),
            "faces": int(number(row, "faces")),
            "watertight": int(number(row, "watertight")),
            "manifold": int(number(row, "manifold")),
            "two_sided_wrap": int(number(row, "two_sided_wrap")),
            "valid": int(number(row, "valid")),
        })
    write_rows(RESULTS / "advanced_alpha_sweep_summary.csv", alpha_rows)

    neighborhood = [entry["volume_ml"] for entry in alpha_rows
                    if entry["alpha_mm"] in (6.25, 7.5) and
                    entry["offset_mm"] in (0.1, 0.2, 0.3)]
    neighborhood_mean = float(np.mean(neighborhood))
    neighborhood_cv = 100.0 * float(np.std(neighborhood)) / neighborhood_mean
    legacy = read_csv(RESULTS / "denoising_benchmark.csv")
    baseline_map: dict[tuple[float, float], float] = {}
    for row in legacy:
        alpha, offset = number(row, "alpha_mm"), number(row, "offset_mm")
        if (row["denoise_method"] == "NONE" and number(row, "voxel_size_mm") == 1.5 and
                number(row, "injected_outlier_percent") == 0 and
                alpha in (6.25, 7.5) and offset in (0.1, 0.2, 0.3)):
            baseline_map.setdefault((alpha, offset), number(row, "volume_ml"))
    baseline_neighborhood = np.array(list(baseline_map.values()))
    baseline_neighborhood_cv = (100.0 * float(np.std(baseline_neighborhood)) /
                                float(np.mean(baseline_neighborhood)))

    stress_rows = [row for row in rows if row["injected_outlier_mode"] == "clustered"]
    by_pipeline: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in stress_rows:
        by_pipeline[row["denoise_pipeline"]].append(row)
    stress_detail: list[dict[str, object]] = []
    stress_summary: list[dict[str, object]] = []
    for pipeline, group in by_pipeline.items():
        baseline = last_match(group, lambda item: number(item, "injected_outlier_percent") == 0)
        base_before = int(number(baseline, "points_before"))
        base_after = int(number(baseline, "points_after"))
        base_volume = number(baseline, "volume_ml")
        method_details = []
        for row in group:
            if number(row, "injected_outlier_percent") == 0:
                continue
            injected = int(number(row, "points_before")) - base_before
            output_excess = int(number(row, "points_after")) - base_after
            estimated_removed = int(np.clip(injected - output_excess, 0, injected))
            recall = 100.0 * estimated_removed / injected if injected else 100.0
            detail = {
                "pipeline": pipeline,
                "noise_percent": number(row, "injected_outlier_percent"),
                "noise_distance_mm": number(row, "injected_outlier_distance_mm"),
                "injected_points": injected,
                "estimated_noise_removed": estimated_removed,
                "estimated_recovery_percent": recall,
                "output_point_excess_vs_baseline": output_excess,
                "volume_delta_ml_vs_method_baseline": number(row, "volume_ml") - base_volume,
                "absolute_volume_delta_ml": abs(number(row, "volume_ml") - base_volume),
                "deployment_total_ms": number(row, "deployment_total_ms"),
                "valid": int(number(row, "valid")),
            }
            stress_detail.append(detail)
            method_details.append(detail)
        stress_summary.append({
            "pipeline": pipeline,
            "baseline_remaining_points": base_after,
            "baseline_volume_ml": base_volume,
            "mean_estimated_recovery_percent": float(np.mean(
                [entry["estimated_recovery_percent"] for entry in method_details])),
            "mean_absolute_volume_delta_ml": float(np.mean(
                [entry["absolute_volume_delta_ml"] for entry in method_details])),
            "max_absolute_volume_delta_ml": float(np.max(
                [entry["absolute_volume_delta_ml"] for entry in method_details])),
            "valid_cases": sum(entry["valid"] for entry in method_details),
            "total_cases": len(method_details),
        })
    stress_summary.sort(key=lambda item: item["pipeline"])
    write_rows(RESULTS / "advanced_cluster_stress_detail.csv", stress_detail)
    write_rows(RESULTS / "advanced_cluster_stress_summary.csv", stress_summary)

    repeats = read_csv(REPEAT)
    metrics = ("voxel_ms", "component_ms", "denoise_ms", "alpha_wrap_ms",
               "validation_ms", "volume_ms", "deployment_total_ms", "research_total_ms")
    repeat_summary: list[dict[str, object]] = []
    for metric in metrics:
        values = np.array([number(row, metric) for row in repeats])
        repeat_summary.append({
            "metric": metric, "n": len(values), "mean": float(np.mean(values)),
            "median": float(np.median(values)), "std": float(np.std(values, ddof=1)),
            "min": float(np.min(values)), "max": float(np.max(values)),
            "p95": float(np.percentile(values, 95)),
        })
    volume_values = np.array([number(row, "volume_ml") for row in repeats])
    volume_cv_repeat = 100.0 * float(np.std(volume_values, ddof=1)) / float(np.mean(volume_values))
    write_rows(RESULTS / "advanced_repeat_statistics.csv", repeat_summary)

    selected = (("BEST_QUALITY", 6.25, 0.1), ("BEST_BALANCED", 7.5, 0.2),
                ("BEST_FAST", 10.0, 0.3))
    recommendations: list[dict[str, object]] = []
    for label, alpha, offset in selected:
        row = alpha_map[(alpha, offset)]
        recommendations.append({
            "candidate": label, "voxel_size_mm": 1.5,
            "pipeline": "VOXEL+COMPONENT", "component_radius_mm": 4.0,
            "alpha_mm": alpha, "offset_mm": offset,
            "remaining_points": int(number(row, "points_after")),
            "cleaning_ms_single_run": number(row, "denoise_ms"),
            "alpha_wrap_ms_single_run": number(row, "alpha_wrap_ms"),
            "deployment_total_ms_single_run": number(row, "deployment_total_ms"),
            "volume_ml": number(row, "volume_ml"),
            "vertices": int(number(row, "vertices")), "faces": int(number(row, "faces")),
            "watertight": int(number(row, "watertight")),
            "manifold": int(number(row, "manifold")),
            "two_sided_wrap": int(number(row, "two_sided_wrap")),
            "valid": int(number(row, "valid")),
        })
    write_rows(RESULTS / "advanced_recommended_candidates.csv", recommendations)

    write_rows(RESULTS / "advanced_key_metrics.csv", [{
        "alpha_neighborhood_definition": "alpha=6.25|7.5; offset=0.1|0.2|0.3",
        "alpha_neighborhood_volume_mean_ml": neighborhood_mean,
        "alpha_neighborhood_volume_cv_percent": neighborhood_cv,
        "voxel_only_matching_neighborhood_cv_percent": baseline_neighborhood_cv,
        "repeat_volume_mean_ml": float(np.mean(volume_values)),
        "repeat_volume_cv_percent": volume_cv_repeat,
        "repeat_deployment_p95_ms": next(item["p95"] for item in repeat_summary
                                         if item["metric"] == "deployment_total_ms"),
        "repeat_alpha_wrap_mean_ms": next(item["mean"] for item in repeat_summary
                                          if item["metric"] == "alpha_wrap_ms"),
    }])

    # Figures
    v15 = [row for row in ablation if row["voxel_size_mm"] == 1.5]
    labels = [row["pipeline"].replace("VOXEL+", "").replace("+", "\n") for row in v15]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)
    x = np.arange(len(v15))
    axes[0].bar(x, [row["cleaning_ms"] for row in v15], color="#457b9d")
    axes[0].set_ylabel("Cleaning time (ms)")
    axes[0].set_title("Denoising ablation at voxel 1.5 mm")
    axes[1].bar(x, [row["volume_ml"] for row in v15], color="#2a9d8f")
    axes[1].set_ylabel("Volume (mL)")
    axes[1].set_title("Same Alpha Wrap: alpha 6.25, offset 0.10 mm")
    for axis in axes:
        axis.set_xticks(x, labels, rotation=0, fontsize=8)
        axis.grid(axis="y", alpha=0.25)
    fig.savefig(FIGURES / "advanced_ablation.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for offset in sorted({row["offset_mm"] for row in alpha_rows}):
        chosen = [row for row in alpha_rows if row["offset_mm"] == offset]
        axes[0].plot([row["alpha_mm"] for row in chosen],
                     [row["volume_ml"] for row in chosen], marker="o", label=f"{offset:g} mm")
        axes[1].plot([row["alpha_mm"] for row in chosen],
                     [row["alpha_wrap_ms"] for row in chosen], marker="o", label=f"{offset:g} mm")
    axes[0].set(title="Alpha/offset volume sensitivity", xlabel="Alpha (mm)", ylabel="Volume (mL)")
    axes[1].set(title="Alpha Wrap runtime", xlabel="Alpha (mm)", ylabel="Runtime (ms)")
    axes[0].legend(title="Offset", ncol=2, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURES / "advanced_alpha_sensitivity.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    for pipeline in sorted(by_pipeline):
        chosen = [row for row in stress_detail if row["pipeline"] == pipeline and
                  row["noise_percent"] == 1.0]
        chosen.sort(key=lambda item: item["noise_distance_mm"])
        axes[0].plot([row["noise_distance_mm"] for row in chosen],
                     [row["estimated_recovery_percent"] for row in chosen], marker="o", label=pipeline)
        axes[1].plot([row["noise_distance_mm"] for row in chosen],
                     [row["absolute_volume_delta_ml"] for row in chosen], marker="o", label=pipeline)
    axes[0].set(title="Cluster-noise recovery (1% injection)", xlabel="Distance (mm)",
                ylabel="Estimated recovery (%)", ylim=(-3, 103))
    axes[1].set(title="Volume restoration (1% injection)", xlabel="Distance (mm)",
                ylabel="Absolute volume change (mL)")
    axes[0].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURES / "advanced_cluster_stress.png", dpi=180)
    plt.close(fig)

    metric_order = ("voxel_ms", "component_ms", "alpha_wrap_ms", "deployment_total_ms")
    fig, axis = plt.subplots(figsize=(9, 5), constrained_layout=True)
    values = [[number(row, metric) for row in repeats] for metric in metric_order]
    axis.boxplot(values, tick_labels=("Voxel", "Component", "Alpha Wrap", "Deployment total"),
                 showmeans=True)
    axis.set_ylabel("Runtime (ms)")
    axis.set_title("Best-overall pipeline, 20 measured runs (warm-up excluded)")
    axis.grid(axis="y", alpha=0.25)
    fig.savefig(FIGURES / "advanced_runtime_repeat.png", dpi=180)
    plt.close(fig)

    print(f"Analyzed {len(rows)} full rows and {len(repeats)} repeat rows")
    print(f"Reasonable-neighborhood volume CV: {neighborhood_cv:.4f}%")
    print(f"20-run volume CV: {volume_cv_repeat:.8f}%")


if __name__ == "__main__":
    main()
