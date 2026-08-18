#!/usr/bin/env python3
"""Measure robust transverse dimensions in 2 mm bins at both cucumber ends."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
POINTS = ROOT / "results" / "denoised_points"
OUTPUT = ROOT / "results" / "end_cross_section_metrics.csv"
METHODS = ("none", "sor", "ror", "sor_ror")
BINS = tuple((float(start), float(start + 2)) for start in range(0, 10, 2))


def read_binary_xyz_ply(path: Path) -> np.ndarray:
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
    match = re.search(r"element vertex (\d+)", header)
    if not match or "format binary_little_endian 1.0" not in header:
        raise RuntimeError(f"Unsupported PLY format: {path}")
    count = int(match.group(1))
    if count == 0:
        return np.empty((0, 3), dtype=float)
    values = np.frombuffer(raw, dtype="<f8", count=count * 3, offset=header_end)
    return values.reshape(count, 3).copy()


def robust_section(points: np.ndarray, center: np.ndarray,
                   transverse: np.ndarray) -> tuple[int, float, float, float]:
    if len(points) < 3:
        return len(points), float("nan"), float("nan"), float("nan")
    uv = (points - center) @ transverse
    low, high = np.percentile(uv, [2.5, 97.5], axis=0)
    width, height = high - low
    return len(points), float(width), float(height), float(np.pi * width * height / 4.0)


def main() -> None:
    baseline = read_binary_xyz_ply(POINTS / "none_remaining.ply")
    center = np.mean(baseline, axis=0)
    covariance = np.cov(baseline - center, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)
    axis = eigenvectors[:, order[-1]]
    if axis[np.argmax(np.abs(axis))] < 0:
        axis = -axis
    transverse = eigenvectors[:, order[:2]]
    baseline_s = (baseline - center) @ axis
    s_min, s_max = float(np.min(baseline_s)), float(np.max(baseline_s))

    outputs: list[dict[str, object]] = []
    for method in METHODS:
        after = read_binary_xyz_ply(POINTS / f"{method}_remaining.ply")
        removed = read_binary_xyz_ply(POINTS / f"{method}_removed.ply")
        datasets = (("before", baseline), ("after", after))
        projections = {stage: (points - center) @ axis for stage, points in datasets}
        removed_s = (removed - center) @ axis if len(removed) else np.empty(0)
        for end in ("bottom", "top"):
            for start, stop in BINS:
                reference = {}
                for stage, points in datasets:
                    s = projections[stage]
                    distance = s - s_min if end == "bottom" else s_max - s
                    mask = (distance >= start) & (distance < stop)
                    reference[stage] = robust_section(points[mask], center, transverse)
                removed_distance = (removed_s - s_min if end == "bottom"
                                    else s_max - removed_s)
                removed_count = int(np.sum((removed_distance >= start) &
                                           (removed_distance < stop)))
                before = reference["before"]
                for stage in ("before", "after"):
                    count, width, height, area = reference[stage]
                    def delta(value: float, base: float) -> float:
                        if stage == "before" or not np.isfinite(base) or base == 0:
                            return 0.0
                        return 100.0 * (value - base) / base
                    outputs.append({
                        "method": method.upper(),
                        "end": end.upper(),
                        "bin_start_mm": start,
                        "bin_end_mm": stop,
                        "stage": stage.upper(),
                        "point_count": count,
                        "cross_section_width_mm": width,
                        "cross_section_height_mm": height,
                        "ellipse_area_mm2": area,
                        "width_change_percent": delta(width, before[1]),
                        "height_change_percent": delta(height, before[2]),
                        "area_change_percent": delta(area, before[3]),
                        "removed_points_in_bin": removed_count if stage == "after" else 0,
                    })

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(outputs[0]))
        writer.writeheader()
        writer.writerows(outputs)
    print(f"Wrote {len(outputs)} rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
