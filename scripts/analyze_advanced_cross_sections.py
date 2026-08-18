#!/usr/bin/env python3
"""Compare residual-noise cleaning stages in 2 mm bins at both ends."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CLEANING = ROOT / "results" / "residual_cleaning"
OUTPUT = ROOT / "results" / "advanced_end_cross_section_metrics.csv"
STAGES = {
    "SOR_ROR": "00_original.ply",
    "COMPONENT": "04_component_remaining.ply",
    "ADAPTIVE_LOCAL": "05_local_remaining.ply",
}
BINS = tuple((float(start), float(start + 2)) for start in range(0, 10, 2))


def read_xyz_ply(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    marker = b"end_header\n"
    offset = raw.find(marker)
    if offset < 0:
        marker = b"end_header\r\n"
        offset = raw.find(marker)
    if offset < 0:
        raise RuntimeError(f"PLY header terminator not found: {path}")
    header_end = offset + len(marker)
    header = raw[:header_end].decode("ascii")
    count_match = re.search(r"element vertex (\d+)", header)
    if not count_match:
        raise RuntimeError(f"Missing vertex count: {path}")
    count = int(count_match.group(1))
    if count == 0:
        return np.empty((0, 3), dtype=float)
    if "format binary_little_endian 1.0" in header:
        # Project point writers store xyz as doubles, optionally followed by RGB.
        vertex_lines = header.split("element vertex", 1)[1].split("element", 1)[0]
        has_rgb = "property uchar red" in vertex_lines
        if has_rgb:
            dtype = np.dtype([("x", "<f8"), ("y", "<f8"), ("z", "<f8"),
                              ("r", "u1"), ("g", "u1"), ("b", "u1")])
            values = np.frombuffer(raw, dtype=dtype, count=count, offset=header_end)
            return np.column_stack((values["x"], values["y"], values["z"]))
        return np.frombuffer(raw, dtype="<f8", count=count * 3,
                             offset=header_end).reshape(count, 3).copy()
    raise RuntimeError(f"Unsupported PLY format: {path}")


def section_metrics(points: np.ndarray, center: np.ndarray,
                    transverse: np.ndarray) -> tuple[int, float, float, float]:
    if len(points) < 3:
        return len(points), float("nan"), float("nan"), float("nan")
    uv = (points - center) @ transverse
    low, high = np.percentile(uv, [2.5, 97.5], axis=0)
    dimensions = np.sort(high - low)[::-1]
    width, height = map(float, dimensions)
    return len(points), width, height, float(np.pi * width * height / 4.0)


def change_percent(value: float, baseline: float) -> float:
    if not np.isfinite(value) or not np.isfinite(baseline) or baseline == 0:
        return float("nan")
    return 100.0 * (value - baseline) / baseline


def main() -> None:
    clouds = {name: read_xyz_ply(CLEANING / filename)
              for name, filename in STAGES.items()}
    baseline = clouds["SOR_ROR"]
    center = np.mean(baseline, axis=0)
    covariance = np.cov(baseline - center, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)
    axis = eigenvectors[:, order[-1]]
    if axis[np.argmax(np.abs(axis))] < 0:
        axis = -axis
    transverse = eigenvectors[:, order[:2]]
    # Use the visually confirmed main component endpoints as the section
    # origin.  The noisy baseline extrema are themselves the disconnected
    # clusters under investigation and would make the first bins describe
    # only noise rather than cucumber geometry.
    main_projection = (clouds["COMPONENT"] - center) @ axis
    s_min, s_max = float(np.min(main_projection)), float(np.max(main_projection))

    rows: list[dict[str, object]] = []
    for end in ("BOTTOM", "TOP"):
        for start, stop in BINS:
            metrics: dict[str, tuple[int, float, float, float]] = {}
            for name, points in clouds.items():
                projection = (points - center) @ axis
                distance = projection - s_min if end == "BOTTOM" else s_max - projection
                selected = points[(distance >= start) & (distance < stop)]
                metrics[name] = section_metrics(selected, center, transverse)
            before = metrics["SOR_ROR"]
            for name in STAGES:
                count, width, height, area = metrics[name]
                rows.append({
                    "endpoint_reference": "CONFIRMED_MAIN_COMPONENT",
                    "stage": name,
                    "end": end,
                    "bin_start_mm": start,
                    "bin_end_mm": stop,
                    "point_count": count,
                    "removed_vs_before": before[0] - count,
                    "width_mm": width,
                    "height_mm": height,
                    "ellipse_area_mm2": area,
                    "width_change_percent": change_percent(width, before[1]),
                    "height_change_percent": change_percent(height, before[2]),
                    "area_change_percent": change_percent(area, before[3]),
                })

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
