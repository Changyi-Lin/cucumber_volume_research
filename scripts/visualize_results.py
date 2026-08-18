#!/usr/bin/env python3
"""Create point-cloud, mesh, overlay, end-cap, and PCA cross-section figures."""

from __future__ import annotations

from pathlib import Path
import struct
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "results" / "figures"
MESH_PATH = ROOT / "results" / "meshes" / "best_balanced_mesh.ply"
RAW_PATH = ROOT.parent / "pointcloud_3100.ply"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from inspect_mesh_components import read_cgal_ply  # noqa: E402


PLY_DTYPES = {
    "char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1",
    "short": "<i2", "int16": "<i2", "ushort": "<u2", "uint16": "<u2",
    "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4",
    "float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8",
}


def read_raw_ply(path: Path) -> tuple[np.ndarray, np.ndarray]:
    properties = []
    vertex_count = 0
    current_element = None
    with path.open("rb") as handle:
        while True:
            line = handle.readline().decode("ascii").strip()
            parts = line.split()
            if parts[:1] == ["format"] and parts[1] != "binary_little_endian":
                raise ValueError("Visualization reader expects binary_little_endian PLY")
            if parts[:1] == ["element"]:
                current_element = parts[1]
                if current_element == "vertex":
                    vertex_count = int(parts[2])
            if parts[:1] == ["property"] and current_element == "vertex":
                properties.append((parts[2], PLY_DTYPES[parts[1]]))
            if line == "end_header":
                data_offset = handle.tell()
                break
    data = np.memmap(path, dtype=np.dtype(properties), mode="r", offset=data_offset,
                     shape=(vertex_count,))
    xyz = np.column_stack((data["x"], data["y"], data["z"])).astype(np.float64)
    if all(name in data.dtype.names for name in ("red", "green", "blue")):
        rgb = np.column_stack((data["red"], data["green"], data["blue"])).astype(np.float32) / 255.0
    else:
        rgb = np.full((vertex_count, 3), 0.25, dtype=np.float32)
    return xyz, rgb


def pca_frame(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = points.mean(axis=0)
    covariance = np.cov(points - mean, rowvar=False, bias=True)
    _, axes = np.linalg.eigh(covariance)
    return mean, axes


def transform(points: np.ndarray, mean: np.ndarray, axes: np.ndarray) -> np.ndarray:
    q = (points - mean) @ axes
    return q[:, [2, 1, 0]]  # length, wide transverse, narrow transverse


def equal_3d(ax, q: np.ndarray) -> None:
    span = np.ptp(q, axis=0)
    ax.set_box_aspect(span)
    ax.set_xlabel("Principal axis (mm)")
    ax.set_ylabel("Transverse 1 (mm)")
    ax.set_zlabel("Transverse 2 (mm)")


def save(fig, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plane_segments(vertices: np.ndarray, faces: np.ndarray, x_value: float) -> list[np.ndarray]:
    segments = []
    for triangle in vertices[faces]:
        hits = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = triangle[i], triangle[j]
            da, db = a[0] - x_value, b[0] - x_value
            if da * db < 0 or abs(da) < 1e-12 or abs(db) < 1e-12:
                if abs(b[0] - a[0]) < 1e-12:
                    continue
                t = (x_value - a[0]) / (b[0] - a[0])
                if -1e-10 <= t <= 1 + 1e-10:
                    hit = a + t * (b - a)
                    if not any(np.linalg.norm(hit - old) < 1e-8 for old in hits):
                        hits.append(hit)
        if len(hits) == 2:
            segments.append(np.asarray(hits))
    return segments


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    raw, rgb = read_raw_ply(RAW_PATH)
    mesh_vertices, mesh_faces = read_cgal_ply(MESH_PATH)
    indices = np.linspace(0, len(raw) - 1, min(100000, len(raw)), dtype=int)
    raw_sample, color_sample = raw[indices], rgb[indices]
    mean, axes = pca_frame(raw_sample)
    raw_q = transform(raw_sample, mean, axes)
    mesh_q = transform(mesh_vertices, mean, axes)
    display_indices = np.linspace(0, len(raw_q) - 1, min(45000, len(raw_q)), dtype=int)

    fig = plt.figure(figsize=(9.0, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*raw_q[display_indices].T, c=color_sample[display_indices], s=0.35, alpha=0.65,
               rasterized=True)
    equal_3d(ax, raw_q)
    ax.set_title("Original Cucumber Point Cloud (PCA coordinates)")
    ax.view_init(elev=22, azim=-63)
    save(fig, "raw_point_cloud.png")

    fig = plt.figure(figsize=(9.0, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    collection = Poly3DCollection(mesh_q[mesh_faces], facecolor="#67a9cf", edgecolor="#24576d",
                                  linewidth=0.16, alpha=0.95)
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(mesh_q[:, 0], mesh_q[:, 1], mesh_q[:, 2])
    equal_3d(ax, mesh_q)
    ax.set_title("Best Balanced CGAL Alpha Wrap Mesh")
    ax.view_init(elev=22, azim=-63)
    save(fig, "best_balanced_alpha_wrap_mesh.png")

    fig = plt.figure(figsize=(9.0, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*raw_q[display_indices].T, c=color_sample[display_indices], s=0.25, alpha=0.35,
               rasterized=True)
    collection = Poly3DCollection(mesh_q[mesh_faces], facecolor="#fdae61", edgecolor="#a14b20",
                                  linewidth=0.12, alpha=0.30)
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(raw_q[:, 0], raw_q[:, 1], raw_q[:, 2])
    equal_3d(ax, raw_q)
    ax.set_title("Point Cloud + Best Balanced Alpha Wrap Overlay")
    ax.view_init(elev=22, azim=-63)
    save(fig, "balanced_point_mesh_overlay.png")

    fig, axes_plot = plt.subplots(1, 2, figsize=(9.6, 4.6))
    low, high = np.min(raw_q[:, 0]), np.max(raw_q[:, 0])
    for ax, side in zip(axes_plot, ("Bottom", "Top")):
        if side == "Bottom":
            tip_mask = raw_q[:, 0] < low + 5
            inner_mask = (raw_q[:, 0] > low + 15) & (raw_q[:, 0] < low + 25)
            end_mask = raw_q[:, 0] < low + 32
        else:
            tip_mask = raw_q[:, 0] > high - 5
            inner_mask = (raw_q[:, 0] < high - 15) & (raw_q[:, 0] > high - 25)
            end_mask = raw_q[:, 0] > high - 32
        local_axis = raw_q[tip_mask].mean(axis=0) - raw_q[inner_mask].mean(axis=0)
        local_axis /= np.linalg.norm(local_axis)
        transverse_1 = np.array([0.0, 1.0, 0.0])
        transverse_1 -= transverse_1.dot(local_axis) * local_axis
        transverse_1 /= np.linalg.norm(transverse_1)
        transverse_2 = np.cross(local_axis, transverse_1)
        local_frame = np.column_stack((local_axis, transverse_1, transverse_2))
        raw_local = raw_q @ local_frame
        mesh_local = mesh_q @ local_frame
        x_value = np.max(raw_local[end_mask, 0]) - 3.0
        raw_mask = end_mask & (np.abs(raw_local[:, 0] - x_value) < 0.8)
        ax.scatter(raw_local[raw_mask, 1], raw_local[raw_mask, 2], s=1.0, c="#6a6a6a", alpha=0.35,
                   label="Point cloud")
        for segment in plane_segments(mesh_local, mesh_faces, x_value):
            ax.plot(segment[:, 1], segment[:, 2], color="#d95f02", linewidth=1.1)
        ax.set_title(f"{side} closure (local plane, 3 mm from end)")
        ax.set_xlabel("Transverse 1 (mm)"); ax.set_ylabel("Transverse 2 (mm)")
        ax.set_aspect("equal", adjustable="box"); ax.grid(alpha=0.2)
    save(fig, "top_bottom_closure_zoom.png")

    fractions = [0.20, 0.35, 0.50, 0.65, 0.80]
    fig, axes_plot = plt.subplots(1, 5, figsize=(15.0, 3.2), sharex=True, sharey=True)
    x_min, x_max = np.min(raw_q[:, 0]), np.max(raw_q[:, 0])
    for ax, fraction in zip(axes_plot, fractions):
        x_value = x_min + fraction * (x_max - x_min)
        slab = np.abs(raw_q[:, 0] - x_value) < 0.8
        ax.scatter(raw_q[slab, 1], raw_q[slab, 2], s=1.0, color="#888888", alpha=0.28,
                   rasterized=True)
        segments = plane_segments(mesh_q, mesh_faces, x_value)
        for segment in segments:
            ax.plot(segment[:, 1], segment[:, 2], color="#d95f02", linewidth=1.0)
        ax.set_title(f"{fraction:.0%} length")
        ax.set_aspect("equal", adjustable="box"); ax.grid(alpha=0.18)
        ax.set_xlabel("T1 (mm)")
    axes_plot[0].set_ylabel("T2 (mm)")
    fig.suptitle("PCA Cross-sections: Point Slab and Alpha Wrap Contour", y=1.02)
    save(fig, "balanced_cross_sections.png")


if __name__ == "__main__":
    main()
