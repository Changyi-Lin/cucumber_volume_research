#!/usr/bin/env python3
"""Generate publication-ready denoising benchmark figures."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_end_cross_sections import read_binary_xyz_ply


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
POINTS = RESULTS / "denoised_points"
ORDER = ("NONE", "SOR", "ROR", "SOR_ROR")
LABELS = {"NONE": "Voxel only", "SOR": "SOR", "ROR": "ROR", "SOR_ROR": "SOR + ROR"}
COLORS = {"NONE": "#6b7280", "SOR": "#2563eb", "ROR": "#059669", "SOR_ROR": "#d97706"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def values(rows: list[dict[str, str]], key: str) -> np.ndarray:
    mapping = {row["method"]: float(row[key]) for row in rows}
    return np.asarray([mapping[method] for method in ORDER])


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def finish(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / name, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def bar_figure(rows: list[dict[str, str]], key: str, title: str, ylabel: str,
               filename: str, error_key: str | None = None,
               limit_500: bool = False, zero: bool = True) -> None:
    data = values(rows, key)
    errors = values(rows, error_key) if error_key else None
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = np.arange(len(ORDER))
    bars = ax.bar(x, data, yerr=errors, capsize=4,
                  color=[COLORS[m] for m in ORDER], width=0.68)
    ax.set_xticks(x, [LABELS[m] for m in ORDER])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if limit_500:
        ax.axhline(500, color="#dc2626", linestyle="--", linewidth=1.5,
                   label="500 ms limit")
        ax.legend(frameon=False)
    if not zero:
        spread = max(float(np.ptp(data)), 0.5)
        ax.set_ylim(float(np.min(data) - 0.25 * spread),
                    float(np.max(data) + 0.35 * spread))
    for bar, value in zip(bars, data):
        ax.annotate(f"{value:.2f}",
                    (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 5), textcoords="offset points", ha="center", va="bottom")
    fig.tight_layout()
    finish(fig, filename)


def method_figures() -> None:
    rows = read_csv(RESULTS / "denoising_method_summary.csv")
    bar_figure(rows, "denoise_ms_mean", "Denoising runtime (10-run mean)", "Runtime (ms)",
               "denoising_runtime.png", "denoise_ms_std")
    bar_figure(rows, "removed_percent", "Points removed at voxel = 1.5 mm",
               "Removed points (%)", "points_removed.png")
    bar_figure(rows, "volume_ml_fixed_mean", "Volume at alpha = 6.25 mm, offset = 0.1 mm",
               "Volume (mL)", "volume_by_denoising_method.png", zero=False)
    bar_figure(rows, "alpha_wrap_ms_mean", "Alpha Wrap runtime (10-run mean)",
               "Runtime (ms)", "alpha_wrap_runtime_by_denoising.png", "alpha_wrap_ms_std")
    bar_figure(rows, "volume_cv_mean_within_voxel_percent",
               "Volume stability across alpha / offset neighborhoods",
               "Mean within-voxel CV (%)", "volume_cv_by_denoising.png")
    bar_figure(rows, "deployment_total_ms_mean", "Deployment pipeline runtime (10-run mean)",
               "Runtime (ms)", "total_runtime_by_denoising.png",
               "deployment_total_ms_std", limit_500=True)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4))
    specs = [
        ("denoise_ms_mean", "Denoise runtime", "ms"),
        ("removed_percent", "Points removed", "%"),
        ("volume_ml_fixed_mean", "Volume", "mL"),
        ("alpha_wrap_ms_mean", "Alpha Wrap runtime", "ms"),
        ("volume_cv_mean_within_voxel_percent", "Volume CV", "%"),
        ("deployment_total_ms_mean", "Deployment total", "ms"),
    ]
    x = np.arange(len(ORDER))
    for ax, (key, title, unit) in zip(axes.flat, specs):
        data = values(rows, key)
        ax.bar(x, data, color=[COLORS[m] for m in ORDER], width=0.68)
        ax.set_xticks(x, [LABELS[m] for m in ORDER], rotation=15, ha="right")
        ax.set_title(title)
        ax.set_ylabel(unit)
        if key == "deployment_total_ms_mean":
            ax.axhline(500, color="#dc2626", linestyle="--", linewidth=1.3)
        for i, value in enumerate(data):
            ax.annotate(f"{value:.2f}", (i, value), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8)
    fig.suptitle("Denoising and Alpha Wrap benchmark summary", fontsize=14)
    fig.tight_layout()
    finish(fig, "denoising_method_dashboard.png")


def stress_figure() -> None:
    rows = read_csv(RESULTS / "outlier_stress_analysis.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    for ax, distance in zip(axes, (5.0, 10.0, 20.0)):
        for method in ORDER:
            subset = sorted((row for row in rows
                             if row["method"] == method and
                             float(row["outlier_distance_mm"]) == distance),
                            key=lambda row: float(row["outlier_percent"]))
            x = [0.0] + [float(row["outlier_percent"]) for row in subset]
            y = [0.0] + [float(row["volume_delta_percent"]) for row in subset]
            ax.plot(x, y, marker="o", linewidth=1.8, markersize=4,
                    color=COLORS[method], label=LABELS[method])
        ax.axhline(0, color="#374151", linewidth=0.8)
        ax.set_title(f"Outliers {distance:g} mm outside AABB")
        ax.set_xlabel("Injected outliers (%)")
    axes[0].set_ylabel("Volume change from clean case (%)")
    axes[-1].legend(frameon=False, loc="best")
    fig.suptitle("Artificial outlier stress test")
    fig.tight_layout()
    finish(fig, "outlier_stress_volume.png")


def density_figure() -> None:
    rows = read_csv(RESULTS / "density_analysis.csv")
    regions = ("BOTTOM", "MIDDLE", "TOP")
    x = np.arange(len(ORDER))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.3, 4.6))
    region_colors = ("#7c3aed", "#0891b2", "#db2777")
    for index, (region, color) in enumerate(zip(regions, region_colors)):
        data = []
        for method in ORDER:
            row = next(row for row in rows
                       if row["denoise_method"] == method and row["region"] == region)
            data.append(float(row["median"]))
        ax.bar(x + (index - 1) * width, data, width, label=region.title(), color=color)
    ax.set_xticks(x, [LABELS[m] for m in ORDER])
    ax.set_ylabel("Median nearest-neighbor distance (mm)")
    ax.set_title("Regional point density after preprocessing")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    finish(fig, "regional_density.png")


def end_protection_figure() -> None:
    rows = read_csv(RESULTS / "end_protection_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    for method in ("SOR", "ROR", "SOR_ROR"):
        subset = sorted((row for row in rows if row["denoise_method"] == method),
                        key=lambda row: float(row["end_protection"]))
        x = [100 * float(row["end_protection"]) for row in subset]
        axes[0].plot(x, [float(row["length_change_mm"]) for row in subset],
                     marker="o", color=COLORS[method], label=LABELS[method])
        axes[1].plot(x, [float(row["removed_percent"]) for row in subset],
                     marker="o", color=COLORS[method], label=LABELS[method])
    axes[0].axhline(1.0, color="#dc2626", linestyle="--", linewidth=1.2,
                    label="1 mm warning")
    axes[0].set_ylabel("Length loss (mm)")
    axes[1].set_ylabel("Points removed (%)")
    for ax in axes:
        ax.set_xlabel("Protected length at each end (%)")
        ax.set_xticks([0, 5, 7.5, 10])
    axes[0].set_title("Endpoint length preservation")
    axes[1].set_title("Filtering strength")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    finish(fig, "end_protection_effect.png")


def pca_projection(points: np.ndarray, center: np.ndarray,
                   basis: np.ndarray) -> np.ndarray:
    return (points - center) @ basis


def sample(points: np.ndarray, maximum: int = 12000) -> np.ndarray:
    if len(points) <= maximum:
        return points
    return points[::math.ceil(len(points) / maximum)]


def point_cloud_figures() -> None:
    baseline = read_binary_xyz_ply(POINTS / "none_remaining.ply")
    center = np.mean(baseline, axis=0)
    eigenvalues, eigenvectors = np.linalg.eigh(np.cov(baseline - center, rowvar=False))
    order = np.argsort(eigenvalues)
    axis = eigenvectors[:, order[-1]]
    radial = eigenvectors[:, order[-2]]
    basis = np.column_stack((axis, radial))
    if axis[np.argmax(np.abs(axis))] < 0:
        basis[:, 0] *= -1
    base_xy = pca_projection(baseline, center, basis)

    fig, axes = plt.subplots(4, 3, figsize=(12, 13), sharex=True, sharey=True)
    for row_index, method in enumerate(ORDER):
        remaining = read_binary_xyz_ply(POINTS / f"{method.lower()}_remaining.ply")
        removed = read_binary_xyz_ply(POINTS / f"{method.lower()}_removed.ply")
        remaining_xy = pca_projection(remaining, center, basis)
        removed_xy = pca_projection(removed, center, basis) if len(removed) else np.empty((0, 2))
        datasets = (sample(base_xy), sample(remaining_xy), removed_xy)
        titles = ("Downsampled input", "Remaining", "Removed")
        colors = ("#6b7280", COLORS[method], "#dc2626")
        for col, (ax, data, title, color) in enumerate(zip(axes[row_index], datasets, titles, colors)):
            if col == 2:
                ax.scatter(sample(base_xy, 7000)[:, 0], sample(base_xy, 7000)[:, 1],
                           s=0.4, color="#d1d5db", alpha=0.35, rasterized=True)
            if len(data):
                ax.scatter(data[:, 0], data[:, 1], s=2.0 if col == 2 else 0.5,
                           color=color, alpha=0.8, rasterized=True)
            elif col == 2:
                ax.text(0.5, 0.5, "No points removed", transform=ax.transAxes,
                        ha="center", va="center")
            if row_index == 0:
                ax.set_title(title)
            if col == 0:
                ax.set_ylabel(f"{LABELS[method]}\nRadial axis (mm)")
            if row_index == len(ORDER) - 1:
                ax.set_xlabel("PCA longitudinal axis (mm)")
            ax.set_aspect("equal", adjustable="box")
    fig.suptitle("Original / remaining / removed point projections")
    fig.tight_layout()
    finish(fig, "original_remaining_removed.png")

    s_min, s_max = np.min(base_xy[:, 0]), np.max(base_xy[:, 0])
    fig, axes = plt.subplots(4, 2, figsize=(10, 11.5), sharey=True)
    for row_index, method in enumerate(ORDER):
        remaining = read_binary_xyz_ply(POINTS / f"{method.lower()}_remaining.ply")
        removed = read_binary_xyz_ply(POINTS / f"{method.lower()}_removed.ply")
        remaining_xy = pca_projection(remaining, center, basis)
        removed_xy = pca_projection(removed, center, basis) if len(removed) else np.empty((0, 2))
        for col, (end, lo, hi) in enumerate((("Bottom", s_min, s_min + 20),
                                              ("Top", s_max - 20, s_max))):
            ax = axes[row_index, col]
            base_mask = (base_xy[:, 0] >= lo) & (base_xy[:, 0] <= hi)
            keep_mask = (remaining_xy[:, 0] >= lo) & (remaining_xy[:, 0] <= hi)
            ax.scatter(base_xy[base_mask, 0], base_xy[base_mask, 1], s=1.2,
                       color="#d1d5db", alpha=0.5, rasterized=True)
            ax.scatter(remaining_xy[keep_mask, 0], remaining_xy[keep_mask, 1], s=1.0,
                       color=COLORS[method], alpha=0.75, rasterized=True)
            if len(removed_xy):
                remove_mask = (removed_xy[:, 0] >= lo) & (removed_xy[:, 0] <= hi)
                ax.scatter(removed_xy[remove_mask, 0], removed_xy[remove_mask, 1], s=16,
                           color="#dc2626", marker="x", linewidths=1.0, rasterized=True)
            if row_index == 0:
                ax.set_title(f"{end} 20 mm zoom")
            if col == 0:
                ax.set_ylabel(f"{LABELS[method]}\nRadial axis (mm)")
            if row_index == len(ORDER) - 1:
                ax.set_xlabel("PCA longitudinal axis (mm)")
            ax.set_xlim(lo, hi)
    fig.suptitle("Endpoint zoom: remaining points and removed points (red ×)")
    fig.tight_layout()
    finish(fig, "end_zoom_removed_points.png")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    style()
    method_figures()
    stress_figure()
    density_figure()
    end_protection_figure()
    point_cloud_figures()
    print(f"Figures written to {FIGURES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
