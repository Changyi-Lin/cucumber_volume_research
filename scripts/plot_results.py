#!/usr/bin/env python3
"""Generate the required benchmark and Pareto figures."""

from __future__ import annotations

import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def unique(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    values = {}
    for row in rows:
        key = (round(f(row, "voxel_size_mm"), 6), round(f(row, "alpha_mm"), 6),
               round(f(row, "offset_mm"), 6))
        values[key] = row
    return list(values.values())


def finish(fig: plt.Figure, filename: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def slice_rows(rows, *, alpha=None, offset=None, voxel_min=0.0):
    output = []
    for row in unique(rows):
        if f(row, "voxel_size_mm") < voxel_min:
            continue
        if alpha is not None and abs(f(row, "alpha_mm") - alpha) > 1e-8:
            continue
        if offset is not None and abs(f(row, "offset_mm") - offset) > 1e-8:
            continue
        output.append(row)
    return output


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = load(RESULTS / "alpha_wrap_benchmark.csv")
    points = sorted(load(RESULTS / "point_count_runtime.csv"), key=lambda r: f(r, "input_points"))
    valid = [r for r in unique(rows) if r["valid"] == "1"]

    voxel_rows = sorted(slice_rows(rows, alpha=15.0, offset=0.3, voxel_min=0.5),
                        key=lambda r: f(r, "voxel_size_mm"))
    x = [f(r, "voxel_size_mm") for r in voxel_rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(x, [f(r, "online_core_ms") for r in voxel_rows], "o-", label="Online core")
    ax.plot(x, [f(r, "alpha_wrap_ms") for r in voxel_rows], "s--", label="Alpha Wrap")
    ax.plot(x, [f(r, "downsample_ms") for r in voxel_rows], "^--", label="Voxel downsample")
    ax.axhline(500, color="crimson", linestyle=":", label="500 ms limit")
    ax.set(xlabel="Voxel size (mm)", ylabel="Runtime (ms)", title="Runtime vs Voxel Size (alpha 15 mm, offset 0.3 mm)")
    ax.grid(alpha=0.25); ax.legend()
    finish(fig, "runtime_vs_voxel_size.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(x, [f(r, "volume_ml") for r in voxel_rows], "o-")
    ax.set(xlabel="Voxel size (mm)", ylabel="Volume (mL)", title="Volume vs Voxel Size (alpha 15 mm, offset 0.3 mm)")
    ax.grid(alpha=0.25)
    finish(fig, "volume_vs_voxel_size.png")

    alpha_rows = sorted([r for r in unique(rows)
                         if abs(f(r, "voxel_size_mm") - 2.0) < 1e-8
                         and abs(f(r, "offset_mm") - 0.3) < 1e-8],
                        key=lambda r: f(r, "alpha_mm"))
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    good = [r for r in alpha_rows if r["valid"] == "1"]
    bad = [r for r in alpha_rows if r["valid"] != "1"]
    ax.plot([f(r, "alpha_mm") for r in good], [f(r, "volume_ml") for r in good], "o-", label="Valid")
    if bad:
        ax.scatter([f(r, "alpha_mm") for r in bad], [f(r, "volume_ml") for r in bad],
                   marker="x", s=70, color="crimson", label="Invalid / two-sided")
    ax.set(xlabel="Alpha (mm)", ylabel="Volume (mL)", title="Volume vs Alpha (voxel 2.0 mm, offset 0.3 mm)")
    ax.grid(alpha=0.25); ax.legend()
    finish(fig, "volume_vs_alpha.png")

    offset_rows = sorted([r for r in unique(rows)
                          if abs(f(r, "voxel_size_mm") - 1.5) < 1e-8
                          and abs(f(r, "alpha_mm") - 15.0) < 1e-8],
                         key=lambda r: f(r, "offset_mm"))
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot([f(r, "offset_mm") for r in offset_rows],
            [f(r, "volume_ml") for r in offset_rows], "o-")
    ax.set(xlabel="Offset (mm)", ylabel="Volume (mL)", title="Volume vs Offset (voxel 1.5 mm, alpha 15 mm)")
    ax.grid(alpha=0.25)
    finish(fig, "volume_vs_offset.png")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    scatter = ax.scatter([f(r, "online_core_ms") for r in valid], [f(r, "volume_ml") for r in valid],
                         c=[f(r, "voxel_size_mm") for r in valid], s=24, alpha=0.7, cmap="viridis")
    ax.axvline(500, color="crimson", linestyle=":", label="500 ms limit")
    selected = {"Balanced": (1.0, 6.25, 0.1), "Fastest": (3.0, 40.0, 0.75), "Quality": (0.0, 5.0, 0.1)}
    for label, key in selected.items():
        row = next(r for r in valid if all(abs(f(r, k) - v) < 1e-8
                   for k, v in zip(("voxel_size_mm", "alpha_mm", "offset_mm"), key)))
        ax.scatter(f(row, "online_core_ms"), f(row, "volume_ml"), marker="*", s=150,
                   edgecolor="black", linewidth=0.7, label=label)
    ax.set(xlabel="Online core runtime (ms)", ylabel="Volume (mL)", title="Runtime vs Volume")
    ax.grid(alpha=0.2); ax.legend()
    fig.colorbar(scatter, ax=ax, label="Voxel size (mm)")
    finish(fig, "runtime_vs_volume.png")

    # Local RMS relative volume change over immediate one-parameter neighbors.
    grid = unique([r for r in rows if 0.75 <= f(r, "voxel_size_mm") <= 2.5
                   and 5.0 <= f(r, "alpha_mm") <= 17.5
                   and f(r, "offset_mm") in (0.1, 0.2, 0.3, 0.5)])
    candidates = []
    for row in grid:
        if row["valid"] != "1":
            continue
        key_values = [f(row, k) for k in ("voxel_size_mm", "alpha_mm", "offset_mm")]
        changes = []
        for dimension, field in enumerate(("voxel_size_mm", "alpha_mm", "offset_mm")):
            peers = [other for other in grid if all(
                abs(f(other, ("voxel_size_mm", "alpha_mm", "offset_mm")[j]) - key_values[j]) < 1e-8
                for j in range(3) if j != dimension)]
            peers.sort(key=lambda r: f(r, field))
            index = next(i for i, other in enumerate(peers) if other is row)
            for neighbor_index in (index - 1, index + 1):
                if 0 <= neighbor_index < len(peers):
                    neighbor = peers[neighbor_index]
                    changes.append(100.0 * (f(neighbor, "volume_ml") - f(row, "volume_ml")) /
                                   f(row, "volume_ml"))
        if len(changes) >= 2:
            candidates.append((row, float(np.sqrt(np.mean(np.square(changes))))))

    pareto = []
    for row, stability in candidates:
        runtime = f(row, "online_core_ms")
        if not any(f(other, "online_core_ms") <= runtime and other_stability <= stability
                   and (f(other, "online_core_ms") < runtime or other_stability < stability)
                   for other, other_stability in candidates):
            pareto.append((row, stability))
    pareto.sort(key=lambda item: f(item[0], "online_core_ms"))
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    sc = ax.scatter([f(r, "online_core_ms") for r, s in candidates], [s for r, s in candidates],
                    c=[f(r, "volume_ml") for r, s in candidates], cmap="plasma", s=30, alpha=0.7)
    ax.plot([f(r, "online_core_ms") for r, s in pareto], [s for r, s in pareto],
            "k.--", linewidth=1, label="Pareto front")
    ax.axvline(500, color="crimson", linestyle=":", label="500 ms limit")
    ax.set_yscale("log")
    ax.set(xlabel="Online core runtime (ms)", ylabel="Local volume instability, RMS (%)",
           title="Runtime–Volume Stability Pareto Plot")
    ax.grid(alpha=0.2); ax.legend()
    fig.colorbar(sc, ax=ax, label="Volume (mL)")
    finish(fig, "runtime_volume_stability_pareto.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot([f(r, "input_points") for r in points], [f(r, "alpha_wrap_ms") for r in points], "o-", label="Alpha Wrap")
    ax.plot([f(r, "input_points") for r in points], [f(r, "online_core_ms") for r in points], "s--", label="Online core")
    ax.axhline(500, color="crimson", linestyle=":", label="500 ms limit")
    ax.set(xlabel="Number of input points", ylabel="Runtime (ms)", title="Number of Points vs Alpha Wrap Runtime")
    ax.ticklabel_format(style="sci", axis="x", scilimits=(6, 6))
    ax.grid(alpha=0.25); ax.legend()
    finish(fig, "points_vs_alpha_wrap_runtime.png")


if __name__ == "__main__":
    main()
