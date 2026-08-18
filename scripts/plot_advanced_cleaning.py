#!/usr/bin/env python3
"""Visual QA for residual components and adaptive-local removals."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_end_cross_sections import read_binary_xyz_ply


ROOT = Path(__file__).resolve().parents[1]
CLEANING = ROOT / "results" / "residual_cleaning"
FIGURES = ROOT / "results" / "figures" / "advanced_cleaning"


def read_component_ply(path: Path) -> tuple[np.ndarray, np.ndarray]:
    raw = path.read_bytes()
    marker = b"end_header\n"
    offset = raw.find(marker)
    if offset < 0:
        raise RuntimeError(f"Unsupported component PLY: {path}")
    header_end = offset + len(marker)
    header = raw[:header_end].decode("ascii")
    count = int(re.search(r"element vertex (\d+)", header).group(1))
    dtype = np.dtype([("x", "<f8"), ("y", "<f8"), ("z", "<f8"),
                      ("r", "u1"), ("g", "u1"), ("b", "u1")])
    data = np.frombuffer(raw, dtype=dtype, count=count, offset=header_end)
    xyz = np.column_stack((data["x"], data["y"], data["z"]))
    rgb = np.column_stack((data["r"], data["g"], data["b"])) / 255.0
    return xyz, rgb


def sample(points: np.ndarray, maximum: int = 16000) -> np.ndarray:
    if len(points) <= maximum:
        return points
    return points[::math.ceil(len(points) / maximum)]


def project(points: np.ndarray, center: np.ndarray, basis: np.ndarray) -> np.ndarray:
    return (points - center) @ basis


def finish(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / name, dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    before = read_binary_xyz_ply(CLEANING / "03_ror_remaining.ply")
    component = read_binary_xyz_ply(CLEANING / "04_component_remaining.ply")
    local = read_binary_xyz_ply(CLEANING / "05_local_remaining.ply")
    removed_component = read_binary_xyz_ply(CLEANING / "04_component_removed.ply")
    removed_local = read_binary_xyz_ply(CLEANING / "05_local_removed.ply")
    colored, rgb = read_component_ply(CLEANING / "component_colored.ply")

    center = np.mean(before, axis=0)
    eigenvalues, eigenvectors = np.linalg.eigh(np.cov(before - center, rowvar=False))
    order = np.argsort(eigenvalues)
    axis = eigenvectors[:, order[-1]]
    if axis[np.argmax(np.abs(axis))] < 0:
        axis = -axis
    basis = np.column_stack((axis, eigenvectors[:, order[-2]]))
    before_xy = project(before, center, basis)
    component_xy = project(component, center, basis)
    local_xy = project(local, center, basis)
    removed_component_xy = project(removed_component, center, basis)
    removed_local_xy = project(removed_local, center, basis)
    colored_xy = project(colored, center, basis)
    s_min, s_max = float(np.min(before_xy[:, 0])), float(np.max(before_xy[:, 0]))

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.2,
                         "axes.spines.top": False, "axes.spines.right": False})

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    views = (("All components", s_min, s_max),
             ("Bottom 25 mm", s_min, s_min + 25),
             ("Top 25 mm", s_max - 25, s_max))
    for ax, (title, lo, hi) in zip(axes, views):
        mask = (colored_xy[:, 0] >= lo) & (colored_xy[:, 0] <= hi)
        indices = np.flatnonzero(mask)
        if len(indices) > 18000:
            indices = indices[::math.ceil(len(indices) / 18000)]
        ax.scatter(colored_xy[indices, 0], colored_xy[indices, 1], s=1.2,
                   c=rgb[indices], alpha=0.9, rasterized=True)
        ax.set_title(title)
        ax.set_xlabel("PCA longitudinal axis (mm)")
        ax.set_xlim(lo, hi)
    axes[0].set_ylabel("PCA radial axis (mm)")
    fig.suptitle("Connected-component analysis (main body gray; removed clusters colored)")
    fig.tight_layout()
    finish(fig, "component_analysis.png")

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharex=True, sharey=True)
    stages = (("SOR + ROR\n59,781 points", before_xy, "#6b7280"),
              ("After component\n59,759 points", component_xy, "#2563eb"),
              (f"After adaptive local\n{len(local):,} points", local_xy, "#059669"))
    for ax, (title, data, color) in zip(axes[:3], stages):
        data = sample(data)
        ax.scatter(data[:, 0], data[:, 1], s=0.7, color=color, alpha=0.8, rasterized=True)
        ax.set_title(title)
        ax.set_xlabel("Longitudinal axis (mm)")
    base = sample(before_xy)
    axes[3].scatter(base[:, 0], base[:, 1], s=0.5, color="#d1d5db",
                    alpha=0.35, rasterized=True)
    axes[3].scatter(removed_component_xy[:, 0], removed_component_xy[:, 1], s=18,
                    color="#dc2626", marker="x", linewidths=1.2, label="Component")
    axes[3].scatter(removed_local_xy[:, 0], removed_local_xy[:, 1], s=10,
                    facecolors="none", edgecolors="#ea580c", linewidths=0.9,
                    label="Adaptive local")
    axes[3].set_title("Removed points")
    axes[3].set_xlabel("Longitudinal axis (mm)")
    axes[3].legend(frameon=False, loc="best")
    axes[0].set_ylabel("Radial axis (mm)")
    fig.suptitle("Residual-noise cleaning stages, identical PCA view")
    fig.tight_layout()
    finish(fig, "cleaning_stage_comparison.png")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (title, lo, hi) in zip(axes, views):
        mask = (before_xy[:, 0] >= lo) & (before_xy[:, 0] <= hi)
        base_indices = np.flatnonzero(mask)
        if len(base_indices) > 16000:
            base_indices = base_indices[::math.ceil(len(base_indices) / 16000)]
        ax.scatter(before_xy[base_indices, 0], before_xy[base_indices, 1], s=0.7,
                   color="#9ca3af", alpha=0.35, rasterized=True)
        for data, marker, size in ((removed_component_xy, "x", 20),
                                   (removed_local_xy, "o", 12)):
            removed_mask = (data[:, 0] >= lo) & (data[:, 0] <= hi)
            if marker == "x":
                ax.scatter(data[removed_mask, 0], data[removed_mask, 1], s=size,
                           color="#dc2626", marker=marker, linewidths=1.2)
            else:
                ax.scatter(data[removed_mask, 0], data[removed_mask, 1], s=size,
                           facecolors="none", edgecolors="#dc2626", linewidths=0.9)
        ax.set_title(title)
        ax.set_xlim(lo, hi)
        ax.set_xlabel("PCA longitudinal axis (mm)")
    axes[0].set_ylabel("PCA radial axis (mm)")
    fig.suptitle("All removed points (red) over the SOR + ROR cloud (gray)")
    fig.tight_layout()
    finish(fig, "removed_points_all.png")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (title, lo, hi) in zip(axes, views):
        component_mask = (component_xy[:, 0] >= lo) & (component_xy[:, 0] <= hi)
        indices = np.flatnonzero(component_mask)
        if len(indices) > 16000:
            indices = indices[::math.ceil(len(indices) / 16000)]
        ax.scatter(component_xy[indices, 0], component_xy[indices, 1], s=0.7,
                   color="#9ca3af", alpha=0.35, rasterized=True)
        remove_mask = (removed_local_xy[:, 0] >= lo) & (removed_local_xy[:, 0] <= hi)
        ax.scatter(removed_local_xy[remove_mask, 0], removed_local_xy[remove_mask, 1],
                   s=14, facecolors="none", edgecolors="#dc2626", linewidths=0.9)
        ax.set_title(title)
        ax.set_xlim(lo, hi)
        ax.set_xlabel("PCA longitudinal axis (mm)")
    axes[0].set_ylabel("PCA radial axis (mm)")
    fig.suptitle("Adaptive-local outliers (red circles) after component cleaning")
    fig.tight_layout()
    finish(fig, "local_outlier_visualization.png")

    # Final raw-to-voxel pipeline: verify that the additional 111 disconnected
    # voxel representatives are off-body clusters rather than valid surface.
    raw_dir = ROOT / "results" / "cleaning"
    raw_voxel = read_binary_xyz_ply(raw_dir / "01_voxel.ply")
    raw_remaining = read_binary_xyz_ply(raw_dir / "04_component_remaining.ply")
    raw_removed = read_binary_xyz_ply(raw_dir / "04_component_removed.ply")
    raw_center = np.mean(raw_voxel, axis=0)
    _, raw_vectors = np.linalg.eigh(np.cov(raw_voxel - raw_center, rowvar=False))
    raw_axis = raw_vectors[:, -1]
    if raw_axis[np.argmax(np.abs(raw_axis))] < 0:
        raw_axis = -raw_axis
    raw_basis = np.column_stack((raw_axis, raw_vectors[:, -2]))
    raw_xy = project(raw_voxel, raw_center, raw_basis)
    raw_remaining_xy = project(raw_remaining, raw_center, raw_basis)
    raw_removed_xy = project(raw_removed, raw_center, raw_basis)
    raw_min, raw_max = float(np.min(raw_xy[:, 0])), float(np.max(raw_xy[:, 0]))
    raw_views = (("All voxel points", raw_min, raw_max),
                 ("Bottom 25 mm", raw_min, raw_min + 25),
                 ("Top 25 mm", raw_max - 25, raw_max))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (title, lo, hi) in zip(axes, raw_views):
        mask = (raw_xy[:, 0] >= lo) & (raw_xy[:, 0] <= hi)
        indices = np.flatnonzero(mask)
        if len(indices) > 16000:
            indices = indices[::math.ceil(len(indices) / 16000)]
        ax.scatter(raw_xy[indices, 0], raw_xy[indices, 1], s=0.7,
                   color="#9ca3af", alpha=0.35, rasterized=True)
        removed_mask = (raw_removed_xy[:, 0] >= lo) & (raw_removed_xy[:, 0] <= hi)
        ax.scatter(raw_removed_xy[removed_mask, 0], raw_removed_xy[removed_mask, 1],
                   s=17, color="#dc2626", marker="x", linewidths=1.0)
        ax.set(title=title, xlim=(lo, hi), xlabel="PCA longitudinal axis (mm)")
    axes[0].set_ylabel("PCA radial axis (mm)")
    fig.suptitle(f"Final voxel + component cleaning: {len(raw_removed):,} removed, "
                 f"{len(raw_remaining):,} retained")
    fig.tight_layout()
    finish(fig, "final_component_removed.png")
    print(f"Wrote advanced cleaning figures to {FIGURES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
