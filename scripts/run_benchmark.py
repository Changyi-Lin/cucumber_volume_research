#!/usr/bin/env python3
"""Resumable process-isolated benchmark driver with per-case timeouts."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "cucumber_alpha_wrap.exe"
INPUT_ARG = Path("..") / "pointcloud_3100.ply"
RESULTS = ROOT / "results"
LOGS = RESULTS / "logs"
MAIN_CSV = RESULTS / "alpha_wrap_benchmark.csv"
POINT_CSV = RESULTS / "point_count_runtime.csv"


def existing_keys(csv_path: Path) -> set[tuple[float, float, float, int]]:
    if not csv_path.exists():
        return set()
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return {
            (round(float(row["voxel_size_mm"]), 6),
             round(float(row["alpha_mm"]), 6),
             round(float(row["offset_mm"]), 6),
             int(float(row.get("input_points", 0))))
            for row in csv.DictReader(handle)
            if row.get("failure_reason") != "TIMEOUT"
        }


def run_case(*, voxel: float, alpha: float, offset: float, csv_path: Path,
             timeout: float, max_points: int = 0, tag: str = "case") -> str:
    LOGS.mkdir(parents=True, exist_ok=True)
    relative_csv = csv_path.relative_to(ROOT)
    command = [str(EXE), "--input", str(INPUT_ARG), "--voxel", str(voxel),
               "--alpha", str(alpha), "--offset", str(offset),
               "--result-csv", str(relative_csv)]
    if max_points:
        command += ["--max-points", str(max_points)]
    env = os.environ.copy()
    env["PATH"] = str(Path(r"C:\msys64\ucrt64\bin")) + os.pathsep + env.get("PATH", "")
    name = f"{tag}_v{voxel:g}_a{alpha:g}_o{offset:g}_n{max_points or 'all'}"
    started = time.perf_counter()
    try:
        completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   timeout=timeout, check=False)
        output = completed.stdout
        status = "OK" if completed.returncode in (0, 2) else f"EXIT_{completed.returncode}"
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + f"\nTIMEOUT after {timeout:.1f} s\n"
        status = "TIMEOUT"
    elapsed = time.perf_counter() - started
    (LOGS / f"{name}.log").write_text(output, encoding="utf-8", errors="replace")
    print(f"{status:>8} {name:<44} wall={elapsed:8.3f}s", flush=True)
    return status


def phase_coarse(timeout: float, resume: bool) -> None:
    keys = existing_keys(MAIN_CSV) if resume else set()
    for voxel in [1.0, 1.5, 2.0, 2.5, 3.0]:
        # Expected counts from the measured voxel table, used only for resume matching.
        expected = {1.0: 164821, 1.5: 60129, 2.0: 29223, 2.5: 16760, 3.0: 10745}[voxel]
        for alpha in [20.0, 25.0, 30.0, 35.0, 40.0]:
            for offset in [0.2, 0.3, 0.5, 0.75]:
                key = (voxel, alpha, offset, expected)
                if key in keys:
                    print(f"    SKIP coarse v{voxel:g} a{alpha:g} o{offset:g}", flush=True)
                    continue
                run_case(voxel=voxel, alpha=alpha, offset=offset,
                         csv_path=MAIN_CSV, timeout=timeout, tag="coarse")


def phase_offset(timeout: float, resume: bool) -> None:
    keys = existing_keys(MAIN_CSV) if resume else set()
    voxel, alpha, expected = 1.5, 15.0, 60129
    for offset in [0.1, 0.2, 0.3, 0.5, 0.75, 1.0]:
        if (voxel, alpha, offset, expected) in keys:
            print(f"    SKIP offset o{offset:g}", flush=True)
            continue
        run_case(voxel=voxel, alpha=alpha, offset=offset,
                 csv_path=MAIN_CSV, timeout=timeout, tag="offset")


def phase_alpha(timeout: float, resume: bool) -> None:
    keys = existing_keys(MAIN_CSV) if resume else set()
    voxel, offset, expected = 2.0, 0.3, 29223
    for alpha in [15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]:
        if (voxel, alpha, offset, expected) in keys:
            print(f"    SKIP alpha a{alpha:g}", flush=True)
            continue
        run_case(voxel=voxel, alpha=alpha, offset=offset,
                 csv_path=MAIN_CSV, timeout=timeout, tag="alpha")


def phase_voxel(timeout: float, resume: bool) -> None:
    keys = existing_keys(MAIN_CSV) if resume else set()
    counts = {0.5: 736676, 0.75: 322475, 1.0: 164821, 1.25: 95054,
              1.5: 60129, 1.75: 40797, 2.0: 29223, 2.5: 16760, 3.0: 10745}
    for voxel, expected in counts.items():
        key = (voxel, 15.0, 0.3, expected)
        if key in keys:
            print(f"    SKIP voxel v{voxel:g}", flush=True)
            continue
        run_case(voxel=voxel, alpha=15.0, offset=0.3,
                 csv_path=MAIN_CSV, timeout=timeout, tag="voxel")


def phase_fine(timeout: float, resume: bool) -> None:
    keys = existing_keys(MAIN_CSV) if resume else set()
    counts = {1.25: 95054, 1.5: 60129, 1.75: 40797, 2.0: 29223, 2.25: 21884, 2.5: 16760}
    for voxel, expected in counts.items():
        for alpha in [5.0, 6.25, 7.5, 8.75, 10.0, 12.5, 15.0, 17.5]:
            for offset in [0.2, 0.3, 0.5]:
                key = (voxel, alpha, offset, expected)
                if key in keys:
                    print(f"    SKIP fine v{voxel:g} a{alpha:g} o{offset:g}", flush=True)
                    continue
                run_case(voxel=voxel, alpha=alpha, offset=offset,
                         csv_path=MAIN_CSV, timeout=timeout, tag="fine")


def phase_points(timeout: float, resume: bool) -> None:
    done = set()
    if resume and POINT_CSV.exists():
        with POINT_CSV.open(newline="", encoding="utf-8") as handle:
            done = {int(float(row["input_points"])) for row in csv.DictReader(handle)}
    for count in [10000, 20000, 30000, 50000, 75000, 100000, 150000, 300000,
                  500000, 750000, 1000000, 1250000, 1573158]:
        if count in done:
            print(f"    SKIP points n{count}", flush=True)
            continue
        run_case(voxel=0.0, alpha=15.0, offset=0.3, max_points=count,
                 csv_path=POINT_CSV, timeout=timeout, tag="points")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["coarse", "alpha", "offset", "voxel", "fine", "points", "all"],
                        default="all")
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not EXE.exists():
        raise SystemExit(f"Executable not found: {EXE}")
    phases = ["coarse", "alpha", "offset", "voxel", "fine", "points"] if args.phase == "all" else [args.phase]
    for phase in phases:
        print(f"\n=== {phase.upper()} ===", flush=True)
        globals()[f"phase_{phase}"](args.timeout, args.resume)
    return 0


if __name__ == "__main__":
    sys.exit(main())
