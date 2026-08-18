#!/usr/bin/env python3
"""Process-isolated advanced cleaning and Alpha Wrap experiment matrix."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "cucumber_alpha_wrap.exe"
RAW = Path("..") / "pointcloud_3100.ply"
PRECLEANED = Path("results") / "denoised_points" / "sor_ror_remaining.ply"
FULL_CSV = Path("results") / "full_cleaning_alpha_wrap_benchmark.csv"
LOGS = ROOT / "results" / "logs" / "advanced_cleaning"


def name_of(case: dict) -> str:
    fields = [case.get("tag", "case"), case.get("pipeline", case.get("method", "none")),
              f"v{case.get('voxel', 0):g}", f"a{case.get('alpha', 6.25):g}",
              f"o{case.get('offset', .1):g}"]
    for key, prefix in (("cluster_radius", "cr"), ("local_k", "lk"),
                        ("mad", "m"), ("end", "e"), ("end_multiplier", "em"),
                        ("outlier_percent", "p"), ("outlier_distance", "d")):
        if key in case:
            fields.append(f"{prefix}{case[key]:g}")
    return "_".join(fields).replace(".", "p").replace("+", "_")


def command_for(case: dict, result_csv: Path = FULL_CSV) -> list[str]:
    voxel = float(case.get("voxel", 0.0))
    method = case.get("method", "none")
    command = [str(EXE), "--input", str(case.get("input", RAW)),
               "--voxel", str(voxel), "--denoise", method,
               "--sor-k", str(case.get("sor_k", 30)),
               "--sor-std", str(case.get("sor_std", 3.0)),
               "--ror-radius", str(case.get("ror_radius", max(1.0, 2.5 * voxel))),
               "--ror-min", str(case.get("ror_min", 3)),
               "--end-protection", str(case.get("end", 0.10)),
               "--end-region-fraction", str(case.get("metric_end", case.get("end", 0.10))),
               "--alpha", str(case.get("alpha", 6.25)),
               "--offset", str(case.get("offset", 0.1)),
               "--outlier-percent", str(case.get("outlier_percent", 0)),
               "--outlier-distance", str(case.get("outlier_distance", 0)),
               "--outlier-mode", case.get("outlier_mode", "dispersed"),
               "--denoise-result-csv", str(result_csv)]
    if case.get("component"):
        command += ["--clustering", "--cluster-radius", str(case.get("cluster_radius", 4.0)),
                    "--component-strategy", case.get("component_strategy", "largest"),
                    "--component-min-size", str(case.get("component_min_size", 1)),
                    "--component-end-distance-multiplier",
                    str(case.get("component_end_multiplier", 0.0))]
        if case.get("confirmed_component_noise", False):
            command.append("--confirmed-disconnected-noise")
    if case.get("local"):
        command += ["--adaptive-local", "--local-knn", str(case.get("local_k", 20)),
                    "--local-pca-knn", str(case.get("local_pca_k", max(20, case.get("local_k", 20)))),
                    "--local-mad-factor", str(case.get("mad", 5.0)),
                    "--local-end-multiplier", str(case.get("end_multiplier", 1.5)),
                    "--local-score-threshold", str(case.get("score_threshold", 0.8))]
        if not case.get("regional", True):
            command.append("--local-global-threshold")
    if "component_csv" in case:
        command += ["--component-csv", str(case["component_csv"])]
    if "cleaning_dir" in case:
        command += ["--cleaning-output-dir", str(case["cleaning_dir"])]
    if "point_prefix" in case:
        command += ["--output-points-prefix", str(case["point_prefix"])]
    if "mesh_prefix" in case:
        command += ["--output-prefix", str(case["mesh_prefix"])]
    if "density_csv" in case:
        command += ["--density-csv", str(case["density_csv"])]
    if "repeat" in case:
        command += ["--repeat", str(case["repeat"])]
    return command


def run(case: dict, timeout: float = 180.0, result_csv: Path = FULL_CSV) -> str:
    env = os.environ.copy()
    env["PATH"] = str(Path(r"C:\msys64\ucrt64\bin")) + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "8"
    name = name_of(case)
    LOGS.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    try:
        completed = subprocess.run(command_for(case, result_csv), cwd=ROOT, env=env,
                                   text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, timeout=timeout, check=False)
        status = "OK" if completed.returncode in (0, 2) else f"EXIT_{completed.returncode}"
        output = completed.stdout
    except subprocess.TimeoutExpired as exc:
        status = "TIMEOUT"
        output = (exc.stdout or "") + f"\nTIMEOUT after {timeout:.1f}s\n"
    wall = time.perf_counter() - started
    (LOGS / f"{name}.log").write_text(output, encoding="utf-8", errors="replace")
    print(f"{status:>8} {name:<80} {wall:7.3f}s", flush=True)
    return status


def voxel(timeout: float) -> None:
    output = ROOT / "results" / "full_voxel_downsample_benchmark.csv"
    command = [str(EXE), "--input", str(RAW), "--voxel-benchmark-csv",
               str(Path("results") / output.name), "--voxel-list", "1,1.25,1.5,1.75,2,2.5"]
    env = os.environ.copy()
    env["PATH"] = str(Path(r"C:\msys64\ucrt64\bin")) + os.pathsep + env.get("PATH", "")
    completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=timeout, check=False)
    (LOGS / "voxel.log").parent.mkdir(parents=True, exist_ok=True)
    (LOGS / "voxel.log").write_text(completed.stdout, encoding="utf-8", errors="replace")
    print(completed.stdout, end="")


def density(timeout: float) -> None:
    output = ROOT / "results" / "advanced_density_analysis.csv"
    if output.exists():
        output.unlink()
    for label, source in (("raw", RAW), ("sor_ror", PRECLEANED)):
        for end in (0.05, 0.075, 0.10):
            run(dict(tag=f"density_{label}", input=source, method="none", voxel=0,
                     end=0, metric_end=end, density_csv=Path("results") / output.name), timeout)


def component_scan(timeout: float) -> None:
    component_dir = ROOT / "results" / "components"
    component_dir.mkdir(parents=True, exist_ok=True)
    for radius in (2.5, 3.0, 3.5, 4.0, 4.5, 5.0):
        run(dict(tag="component_scan", input=PRECLEANED, method="none", voxel=0,
                 end=0.10, component=True, cluster_radius=radius,
                 component_strategy="largest", component_min_size=1,
                 component_csv=Path("results") / "components" / f"radius_{radius:g}.csv"), timeout)
    for minimum in (1, 3, 5, 10, 20, 50, 100):
        run(dict(tag="component_min", input=PRECLEANED, method="none", voxel=0,
                 end=0.10, component=True, cluster_radius=4.0,
                 component_strategy="min_size", component_min_size=minimum), timeout)


def local_scan(timeout: float) -> None:
    cases = []
    for local_k in (10, 20):
        for mad in (3.0, 4.0, 5.0, 6.0):
            cases.append(dict(local_k=local_k, mad=mad, regional=True,
                              end_multiplier=1.5))
    for local_k in (15, 30):
        for mad in (4.0, 5.0):
            cases.append(dict(local_k=local_k, mad=mad, regional=True,
                              end_multiplier=1.5))
    for mad in (4.0, 5.0):
        cases.append(dict(local_k=20, mad=mad, regional=False, end_multiplier=1.5))
        cases.append(dict(local_k=20, mad=mad, regional=True, end_multiplier=1.25))
    for end_multiplier in (2.0, 3.0):
        cases.append(dict(local_k=10, local_pca_k=20, mad=6.0, regional=True,
                          end_multiplier=end_multiplier))
    for settings in cases:
        run(dict(tag="local_scan", input=PRECLEANED, method="none", voxel=0,
                 end=0.10, component=True, cluster_radius=4.0,
                 component_strategy="largest", local=True, **settings), timeout)


def ablation(timeout: float, local_mad: float = 5.0) -> None:
    for voxel_size in (1.5, 2.0):
        radius = (8.0 / 3.0) * voxel_size
        pipelines = [
            dict(pipeline="A_VOXEL", method="none"),
            dict(pipeline="B_SOR", method="sor"),
            dict(pipeline="C_ROR", method="ror"),
            dict(pipeline="D_SOR_ROR", method="sor_ror"),
            dict(pipeline="E_SOR_ROR_COMPONENT", method="sor_ror", component=True),
            dict(pipeline="F_SOR_ROR_COMPONENT_LOCAL", method="sor_ror", component=True, local=True),
            dict(pipeline="G_COMPONENT_LOCAL", method="none", component=True, local=True),
        ]
        for pipeline in pipelines:
            run(dict(tag="ablation", voxel=voxel_size, end=0.10,
                     ror_radius=2.5 * voxel_size, cluster_radius=radius,
                     confirmed_component_noise=True,
                     local_k=10, local_pca_k=20, mad=local_mad,
                     end_multiplier=3.0,
                     alpha=6.25, offset=0.1, **pipeline), timeout)


def alpha_sweep(timeout: float, local_mad: float = 5.0) -> None:
    # Residual-noise inspection showed that the disconnected components are
    # true noise while adaptive-local candidates lie on the valid surface.
    # Sweep Alpha Wrap on the geometry-preserving component-only candidate.
    base = dict(tag="alpha_sweep", pipeline="VOXEL_COMPONENT", voxel=1.5,
                method="none", component=True, cluster_radius=4.0,
                confirmed_component_noise=True, end=0.10)
    for alpha in (20.0, 25.0, 30.0, 35.0, 40.0, 45.0):
        for offset in (0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
            run(dict(alpha=alpha, offset=offset, **base), timeout)
    for alpha in (6.25, 7.5, 10.0, 15.0, 17.5):
        for offset in (0.1, 0.2, 0.3, 0.5):
            run(dict(tag="alpha_fine", alpha=alpha, offset=offset,
                     **{key: value for key, value in base.items() if key != "tag"}), timeout)


def stress(timeout: float, local_mad: float = 5.0) -> None:
    pipelines = [
        dict(pipeline="SOR_ROR", method="sor_ror"),
        dict(pipeline="COMPONENT", method="none", component=True),
        dict(pipeline="ADAPTIVE_LOCAL", method="none", local=True),
        dict(pipeline="COMPONENT_LOCAL", method="none", component=True, local=True),
    ]
    cases = [(0.0, 0.0)] + [(percent, distance)
                            for percent in (0.01, 0.05, 0.1, 0.5, 1.0)
                            for distance in (3.0, 5.0, 10.0, 20.0)]
    for pipeline in pipelines:
        for percent, distance in cases:
            run(dict(tag="cluster_stress", voxel=1.5, end=0.10,
                     ror_radius=3.75, cluster_radius=4.0,
                     confirmed_component_noise=True,
                     local_k=10, local_pca_k=20, mad=local_mad,
                     end_multiplier=3.0,
                     alpha=6.25, offset=0.1, outlier_percent=percent,
                     outlier_distance=distance, outlier_mode="clustered", **pipeline), timeout)


def outputs(timeout: float, local_mad: float = 5.0) -> None:
    cleaning = Path("results") / "cleaning"
    base = dict(tag="final_output", pipeline="VOXEL_COMPONENT", voxel=1.5,
                method="none", end=0.10, component=True, cluster_radius=4.0,
                confirmed_component_noise=True)
    candidates = (
        ("quality", 6.25, 0.1),
        ("balanced", 7.5, 0.2),
        ("fast", 10.0, 0.3),
    )
    for label, alpha, offset in candidates:
        case = dict(base, tag=f"final_{label}", alpha=alpha, offset=offset,
                    mesh_prefix=Path("results") / "meshes" / f"best_{label}_mesh")
        if label == "quality":
            case.update(cleaning_dir=cleaning,
                        component_csv=cleaning / "components.csv",
                        point_prefix=cleaning / "final")
        run(case, timeout)


def residual_outputs(timeout: float, local_mad: float = 5.0) -> None:
    cleaning = Path("results") / "residual_cleaning"
    run(dict(tag="residual_before", pipeline="SOR_ROR_ONLY", input=PRECLEANED,
             voxel=0, method="none", end=0.10, alpha=6.25, offset=0.1), timeout)
    run(dict(tag="residual_component", pipeline="SOR_ROR_COMPONENT", input=PRECLEANED,
             voxel=0, method="none", end=0.10, component=True, cluster_radius=4.0,
             component_end_multiplier=0.0, confirmed_component_noise=True,
             alpha=6.25, offset=0.1), timeout)
    run(dict(tag="residual_output", pipeline="COMPONENT_LOCAL", input=PRECLEANED,
             voxel=0, method="none", end=0.10, component=True, cluster_radius=4.0,
             confirmed_component_noise=True,
             component_end_multiplier=0.0, local=True, local_k=10, local_pca_k=20,
             mad=local_mad,
             end_multiplier=3.0, alpha=6.25, offset=0.1, cleaning_dir=cleaning,
             component_csv=cleaning / "components.csv",
             point_prefix=cleaning / "final"), timeout)


def repeat(timeout: float, local_mad: float = 5.0) -> None:
    output = Path("results") / "full_cleaning_repeat_best.csv"
    if (ROOT / output).exists():
        (ROOT / output).unlink()
    run(dict(tag="repeat_best", pipeline="VOXEL_COMPONENT", voxel=1.5,
             method="none", end=0.10, component=True, cluster_radius=4.0,
             confirmed_component_noise=True,
             alpha=6.25, offset=0.1, repeat=20), timeout, result_csv=output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("voxel", "density", "component_scan", "local_scan",
                                             "ablation", "alpha_sweep", "stress", "outputs",
                                             "residual_outputs", "repeat", "all"), default="all")
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument("--local-mad", type=float, default=6.0)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()
    if args.fresh and (ROOT / FULL_CSV).exists():
        (ROOT / FULL_CSV).unlink()
    phases = ("voxel", "density", "component_scan", "local_scan", "ablation",
              "alpha_sweep", "stress", "outputs", "residual_outputs", "repeat") \
        if args.phase == "all" else (args.phase,)
    for phase in phases:
        print(f"\n=== {phase.upper()} ===", flush=True)
        function = globals()[phase]
        if phase in ("ablation", "alpha_sweep", "stress", "outputs", "residual_outputs", "repeat"):
            function(args.timeout, args.local_mad)
        else:
            function(args.timeout)


if __name__ == "__main__":
    main()
