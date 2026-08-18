# Cucumber CGAL Alpha Wrap Benchmark

This project reconstructs a watertight cucumber mesh from a PLY point cloud with the real
`CGAL::alpha_wrap_3()` API, validates the mesh, and measures enclosed volume.

## Implemented pipeline

1. Custom ASCII / binary PLY XYZ reader (extra fields are skipped).
2. Packed-integer hash voxel downsampling using voxel centroids.
3. Optional SOR, adaptive ROR, SOR+ROR, endpoint protection, and guarded Euclidean clustering.
4. `CGAL::alpha_wrap_3(point_range, alpha, offset, Surface_mesh)`.
5. Conservative removal of tiny disconnected debris only when the largest raw component
   contains at least 95% of all faces. Raw and cleaned component counts are both recorded.
6. Validation: triangle mesh, closed, connected, combinatorial manifoldness, outward
   orientation, self-intersections, Euler characteristic, genus, area, and PCA line tests for
   an outer+inner two-sided shell.
7. Signed tetrahedral mesh volume; coordinates are interpreted as millimetres and
   `volume_ml = volume_mm3 / 1000`.

`online_core_ms` is the online quantity used for the 500 ms decision:

```text
voxel downsampling + denoising + Alpha Wrap + required component cleanup + volume
```

`reconstruction_only_ms` is exactly `Alpha Wrap + volume`, as requested. `total_ms` also
includes the relatively expensive validation checks. PLY load is reported separately.

## Build (the tested Windows / MSYS2 environment)

Required MSYS2 UCRT64 packages:

```text
mingw-w64-ucrt-x86_64-gcc
mingw-w64-ucrt-x86_64-cmake
mingw-w64-ucrt-x86_64-ninja
mingw-w64-ucrt-x86_64-cgal
mingw-w64-ucrt-x86_64-eigen3
mingw-w64-ucrt-x86_64-tbb
mingw-w64-ucrt-x86_64-python
mingw-w64-ucrt-x86_64-python-numpy
mingw-w64-ucrt-x86_64-python-matplotlib
```

From PowerShell:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
C:\msys64\ucrt64\bin\cmake.exe -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
C:\msys64\ucrt64\bin\cmake.exe --build build --parallel 2
```

The Release build uses `-O3 -march=native -DNDEBUG`.

## Run

Dataset information:

```powershell
.\build\cucumber_alpha_wrap.exe --input ..\pointcloud_3100.ply --info
```

One reconstruction:

```powershell
.\build\cucumber_alpha_wrap.exe --input ..\pointcloud_3100.ply `
  --voxel 1 --alpha 6.25 --offset 0.1 `
  --result-csv results\one_run.csv `
  --output-prefix results\meshes\one_run
```

Repeated benchmark (one warm-up is automatically excluded):

```powershell
.\build\cucumber_alpha_wrap.exe --input ..\pointcloud_3100.ply `
  --voxel 1 --alpha 6.25 --offset 0.1 --repeat 20 `
  --result-csv results\repeat.csv
```

End-protected adaptive ROR and split denoising timings:

```powershell
.\build\cucumber_alpha_wrap.exe --input ..\pointcloud_3100.ply `
  --voxel 1.5 --denoise ror --ror-radius 3.75 --ror-min 3 `
  --end-protection 0.05 --alpha 6.25 --offset 0.1 `
  --denoise-result-csv results\one_denoise_run.csv `
  --output-points-prefix results\denoised_points\ror
```

Filter only (the current research recommendation; no Alpha Wrap):

```powershell
C:\msys64\ucrt64\bin\python.exe filtering.py ..\pointcloud_3100.ply results\cleaned.ply
```

`filtering.py` defaults to `voxel=1.5 mm` and dominant-component filtering at
`4.0 mm` connectivity. It removes disconnected cluster noise only when the
largest component dominates (>=98% of all points, second largest <1%). Use
`--stage-dir results\stages --component-csv results\components.csv` to retain
all intermediate PLYs and component metadata. `--preset local-diagnostic`
enables the conservative research-only local diagnostic (K=10, PCA K=20,
MAD=6); it is not the recommended production default.

Full resumable experiment and plots:

```powershell
C:\msys64\ucrt64\bin\python.exe scripts\run_benchmark.py --phase all --resume --timeout 300
C:\msys64\ucrt64\bin\python.exe scripts\analyze_results.py
C:\msys64\ucrt64\bin\python.exe scripts\plot_results.py
C:\msys64\ucrt64\bin\python.exe scripts\visualize_results.py
C:\msys64\ucrt64\bin\python.exe scripts\run_denoising_benchmark.py --phase all --fresh --timeout 120
C:\msys64\ucrt64\bin\python.exe scripts\analyze_denoising_results.py
C:\msys64\ucrt64\bin\python.exe scripts\analyze_end_cross_sections.py
C:\msys64\ucrt64\bin\python.exe scripts\plot_denoising_results.py
```

## Main outputs

- `results/alpha_wrap_benchmark.csv`: parameter search, topology, volume, and split timing.
- `results/point_count_runtime.csv`: input-point scaling experiment.
- `results/repeat_*.csv` and `repeat_statistics.csv`: repeated-run distributions.
- `results/denoising_benchmark.csv`: SOR/ROR search, endpoint, clustering, and stress runs.
- `results/denoising_repeatability.csv`: 10 interleaved repeats for each preprocessing method.
- `results/denoising_method_summary.csv`, `outlier_stress_*.csv`, and endpoint/density CSVs:
  derived denoising analyses.
- `results/denoised_points/`: representative remaining and removed point PLYs.
- `results/meshes/best_*.ply` and `.stl`: Fastest, Best Balanced, and Best Quality meshes.
- `results/figures/`: required benchmark, Pareto, point-cloud, overlay, end, and cross-section figures.
- `README_RESULTS.md`: the complete measured results and conclusions.
