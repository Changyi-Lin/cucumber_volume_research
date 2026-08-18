#!/usr/bin/env python3
"""Filter a dense cucumber point cloud with the validated production preset.

The default preset implements the current research recommendation:
    voxel 1.5 mm -> dominant connected component (4 mm radius)

It intentionally does not enable SOR, ROR, or Adaptive Local filtering: on the
validated cucumber these either missed cluster noise or removed valid surface
points.  This script delegates the performance-critical work to the project's
Release C++ executable and uses its --filter-only mode, so Alpha Wrap is not
run when a cleaned point cloud is all that is requested.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
DEFAULT_EXE = ROOT / "build" / "cucumber_alpha_wrap.exe"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter a dense cucumber PLY using the validated component-cleaning preset.")
    parser.add_argument("input", type=Path, help="Input binary PLY point cloud.")
    parser.add_argument("output", type=Path, help="Output cleaned PLY point cloud.")
    parser.add_argument("--preset", choices=("recommended", "conservative", "local-diagnostic"),
                        default="recommended", help="Filtering policy (default: recommended).")
    parser.add_argument("--voxel-mm", type=float, default=1.5,
                        help="Voxel size in mm (default: 1.5).")
    parser.add_argument("--component-radius-mm", type=float,
                        help="Connectivity radius; defaults to 2.667 × voxel (4 mm at 1.5 mm).")
    parser.add_argument("--preserve-detached-ends", action="store_true",
                        help="Keep plausible detached end components within 2 × connectivity radius.")
    parser.add_argument("--stage-dir", type=Path,
                        help="Optional directory for all intermediate PLY stages and RGB component PLY.")
    parser.add_argument("--component-csv", type=Path,
                        help="Optional CSV of component ID, size, centroid, AABB, and main-body distance.")
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE,
                        help="Path to cucumber_alpha_wrap.exe.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input PLY not found: {args.input}")
    if args.voxel_mm <= 0:
        raise SystemExit("--voxel-mm must be positive")
    if not args.exe.is_file():
        raise SystemExit(f"Release executable not found: {args.exe}. Build the project first.")

    radius = args.component_radius_mm or (8.0 / 3.0) * args.voxel_mm
    if radius <= 0:
        raise SystemExit("--component-radius-mm must be positive")

    # The dominance guard inside C++ only removes components when the largest
    # component is >=98% of all points and the second largest is <1%.
    end_multiplier = 2.0 if (args.preserve_detached_ends or args.preset == "conservative") else 0.0
    command = [
        str(args.exe), "--input", str(args.input), "--filter-only",
        "--voxel", str(args.voxel_mm), "--denoise", "none",
        "--clustering", "--cluster-radius", str(radius),
        "--component-strategy", "largest", "--component-min-size", "1",
        "--component-end-distance-multiplier", str(end_multiplier),
        "--confirmed-disconnected-noise", "--end-protection", "0.10",
        "--end-region-fraction", "0.10",
        "--output-points-prefix", str(args.output.with_suffix("")),
    ]

    if args.preset == "local-diagnostic":
        command += ["--adaptive-local", "--local-knn", "10", "--local-pca-knn", "20",
                    "--local-mad-factor", "6", "--local-end-multiplier", "3"]
    if args.stage_dir:
        command += ["--cleaning-output-dir", str(args.stage_dir)]
    if args.component_csv:
        command += ["--component-csv", str(args.component_csv)]

    env = os.environ.copy()
    msys_bin = Path(r"C:\msys64\ucrt64\bin")
    if msys_bin.is_dir():
        env["PATH"] = str(msys_bin) + os.pathsep + env.get("PATH", "")
    env.setdefault("OMP_NUM_THREADS", "8")

    completed = subprocess.run(command, cwd=ROOT, env=env, text=True)
    generated = args.output.with_suffix("").with_name(args.output.with_suffix("").name + "_remaining.ply")
    if completed.returncode != 0:
        return completed.returncode
    if not generated.is_file():
        raise SystemExit(f"Filtering completed but expected output was not produced: {generated}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if generated.resolve() != args.output.resolve():
        shutil.move(str(generated), str(args.output))
    removed = generated.with_name(generated.stem.replace("_remaining", "_removed") + ".ply")
    if removed.is_file():
        removed.unlink()
    print(f"Cleaned point cloud: {args.output}")
    print(f"Preset={args.preset}; voxel={args.voxel_mm:g} mm; component radius={radius:g} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
