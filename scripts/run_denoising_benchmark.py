#!/usr/bin/env python3
"""Process-isolated denoising/Alpha-Wrap benchmark matrix."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "cucumber_alpha_wrap.exe"
CSV = Path("results") / "denoising_benchmark.csv"
LOGS = ROOT / "results" / "logs" / "denoising"


def config_name(c: dict) -> str:
    fields = [c.get("tag", "case"), f"v{c['voxel']:g}", c["method"],
              f"k{c.get('k', 30)}", f"s{c.get('std', 2):g}",
              f"r{c.get('radius', 0):g}", f"n{c.get('min_neighbors', 3)}",
              f"e{c.get('end', 0):g}", f"a{c.get('alpha', 6.25):g}",
              f"o{c.get('offset', .1):g}"]
    if c.get("outlier_percent", 0):
        fields += [f"p{c['outlier_percent']:g}", f"d{c['outlier_distance']:g}"]
    if c.get("clustering"):
        fields += ["cluster", f"cr{c['cluster_radius']:g}"]
    return "_".join(fields).replace(".", "p").replace("+", "_")


def run(c: dict, timeout: float = 120.0, extra: list[str] | None = None,
        result_csv: Path = CSV) -> str:
    radius = c.get("radius", 2.5 * c["voxel"])
    command = [str(EXE), "--input", str(Path("..") / "pointcloud_3100.ply"),
               "--voxel", str(c["voxel"]), "--denoise", c["method"],
               "--sor-k", str(c.get("k", 30)), "--sor-std", str(c.get("std", 2.0)),
               "--ror-radius", str(radius), "--ror-min", str(c.get("min_neighbors", 3)),
               "--end-protection", str(c.get("end", 0.0)),
               "--end-region-fraction", str(c.get("metric_end", 0.05)),
               "--alpha", str(c.get("alpha", 6.25)), "--offset", str(c.get("offset", 0.1)),
               "--outlier-percent", str(c.get("outlier_percent", 0)),
               "--outlier-distance", str(c.get("outlier_distance", 0)),
               "--denoise-result-csv", str(result_csv)]
    if c.get("clustering"):
        command += ["--clustering", "--cluster-radius", str(c["cluster_radius"])]
    if extra:
        command += extra
    env = os.environ.copy()
    env["PATH"] = str(Path(r"C:\msys64\ucrt64\bin")) + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "8"
    name = config_name(c)
    LOGS.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    try:
        completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   timeout=timeout, check=False)
        status = "OK" if completed.returncode in (0, 2) else f"EXIT_{completed.returncode}"
        output = completed.stdout
    except subprocess.TimeoutExpired as exc:
        status = "TIMEOUT"
        output = (exc.stdout or "") + f"\nTIMEOUT after {timeout:.1f}s\n"
    wall = time.perf_counter() - started
    (LOGS / f"{name}.log").write_text(output, encoding="utf-8", errors="replace")
    print(f"{status:>8} {name:<72} {wall:7.3f}s", flush=True)
    return status


def coarse(timeout: float) -> None:
    run(dict(tag="coarse", voxel=1.5, method="none"), timeout)
    for end in (0.0, 0.05):
        for k in (20, 30, 50):
            for std in (1.0, 1.5, 2.0, 2.5):
                run(dict(tag="coarse", voxel=1.5, method="sor", k=k, std=std, end=end), timeout)
        for radius in (2.25, 3.0, 3.75, 4.5):
            for minimum in (2, 3, 5, 8):
                run(dict(tag="coarse", voxel=1.5, method="ror", radius=radius,
                         min_neighbors=minimum, end=end), timeout)
    for k, std in ((20, 1.5), (30, 2.0), (30, 2.5), (50, 2.0)):
        for radius, minimum in ((3.75, 3), (4.5, 5)):
            run(dict(tag="coarse", voxel=1.5, method="sor_ror", k=k, std=std,
                     radius=radius, min_neighbors=minimum, end=0.05), timeout)


def fine(timeout: float) -> None:
    """Cover the requested boundary values without an unnecessary full Cartesian grid."""
    for k in (10, 80):
        for std in (0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
            run(dict(tag="fine", voxel=1.5, method="sor", k=k, std=std,
                     end=0.05), timeout)
    for k in (20, 30, 50):
        for std in (0.75, 1.25, 3.0):
            run(dict(tag="fine", voxel=1.5, method="sor", k=k, std=std,
                     end=0.05), timeout)

    ror_cases = set()
    for radius in (2.0, 2.5, 3.0, 4.0):
        for minimum in (2, 3, 5, 8, 10):
            ror_cases.add((radius, minimum))
    ror_cases.update({(1.0, 2), (1.0, 3),
                      (1.5, 2), (1.5, 3), (1.5, 5),
                      (4.0, 15), (5.0, 3), (5.0, 5), (5.0, 8), (5.0, 15)})
    for radius, minimum in sorted(ror_cases):
        run(dict(tag="fine", voxel=1.5, method="ror", radius=radius,
                 min_neighbors=minimum, end=0.05), timeout)


def endpoints(timeout: float) -> None:
    for method in ("sor", "ror", "sor_ror"):
        for end in (0.0, 0.05, 0.075, 0.10):
            run(dict(tag="ends", voxel=1.5, method=method, k=30, std=3.0,
                     radius=3.75, min_neighbors=3, end=end), timeout)


def compare(timeout: float) -> None:
    for voxel in (1.0, 1.5, 2.0, 2.5):
        configurations = [
            dict(method="none", end=0.0),
            dict(method="sor", k=30, std=3.0, end=0.05),
            dict(method="ror", min_neighbors=3, end=0.05),
            dict(method="sor_ror", k=30, std=3.0, min_neighbors=3, end=0.05),
        ]
        for base in configurations:
            for alpha in (6.25, 7.5, 10.0):
                for offset in (0.1, 0.2, 0.3):
                    c = dict(tag="compare", voxel=voxel, radius=2.5 * voxel,
                             alpha=alpha, offset=offset, **base)
                    run(c, timeout)


def clusters(timeout: float) -> None:
    for end in (0.0, 0.05):
        for multiplier in (2.0, 2.5, 3.0):
            run(dict(tag="cluster", voxel=1.5, method="none", end=end,
                     clustering=True, cluster_radius=multiplier * 1.5), timeout)


def stress(timeout: float) -> None:
    configs = [
        dict(method="none", end=0.0),
        dict(method="sor", k=30, std=3.0, end=0.05),
        dict(method="ror", min_neighbors=3, end=0.05),
        dict(method="sor_ror", k=30, std=3.0, min_neighbors=3, end=0.05),
    ]
    cases = [(0.0, 0.0)] + [(percent, distance)
                            for percent in (0.01, 0.05, 0.1, 0.5, 1.0)
                            for distance in (5.0, 10.0, 20.0)]
    for base in configs:
        for percent, distance in cases:
            run(dict(tag="stress", voxel=1.5, radius=3.75, alpha=6.25, offset=0.1,
                     outlier_percent=percent, outlier_distance=distance, **base), timeout)


def outputs(timeout: float) -> None:
    density = Path("results") / "density_analysis.csv"
    if (ROOT / density).exists():
        (ROOT / density).unlink()
    configs = [
        ("none", dict(method="none", end=0.0)),
        ("sor", dict(method="sor", k=30, std=3.0, end=0.05)),
        ("ror", dict(method="ror", min_neighbors=3, end=0.05)),
        ("sor_ror", dict(method="sor_ror", k=30, std=3.0, min_neighbors=3, end=0.05)),
    ]
    for label, base in configs:
        prefix = Path("results") / "denoised_points" / label
        extra = ["--output-points-prefix", str(prefix), "--density-csv", str(density)]
        run(dict(tag="output", voxel=1.5, radius=3.75, alpha=6.25, offset=0.1, **base),
            timeout, extra)


def repeat(timeout: float) -> None:
    repeat_csv = Path("results") / "denoising_repeatability.csv"
    if (ROOT / repeat_csv).exists():
        (ROOT / repeat_csv).unlink()
    configs = [
        dict(method="none", end=0.0),
        dict(method="sor", k=30, std=3.0, end=0.05),
        dict(method="ror", min_neighbors=3, end=0.05),
        dict(method="sor_ror", k=30, std=3.0, min_neighbors=3, end=0.05),
    ]
    for repetition in range(10):
        for base in configs:
            run(dict(tag=f"repeat{repetition + 1:02d}", voxel=1.5, radius=3.75,
                     alpha=6.25, offset=0.1, **base), timeout,
                result_csv=repeat_csv)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["coarse", "fine", "endpoints", "compare", "clusters",
                                                    "stress", "outputs", "repeat", "all"], default="all")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()
    if args.fresh and (ROOT / CSV).exists():
        (ROOT / CSV).unlink()
    phases = ["coarse", "fine", "endpoints", "compare", "clusters", "stress", "outputs", "repeat"] \
        if args.phase == "all" else [args.phase]
    for phase in phases:
        print(f"\n=== {phase.upper()} ===", flush=True)
        globals()[phase](args.timeout)


if __name__ == "__main__":
    main()
