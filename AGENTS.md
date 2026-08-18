# Cucumber volume research handoff

## Research goal

Estimate cucumber volume from a dense PLY point cloud using a geometry-preserving
preprocessing pipeline followed by real `CGAL::alpha_wrap_3`, mesh validation,
and tetrahedral volume integration. The online processing target is under 500 ms
when point data is already in memory.

## Current recommended filtering pipeline

```text
Dense binary PLY
→ voxel downsample (1.5 mm)
→ dominant connected component (radius 4.0 mm)
→ cleaned binary PLY
```

The final production preset intentionally disables SOR, ROR, and Adaptive Local
filtering for the validated cucumber: the residual visible noise is made of
disconnected clusters, while Adaptive Local marked valid surface points.

Use `filtering.py INPUT.ply OUTPUT.ply` for filtering only. It invokes the
Release C++ executable with `--filter-only`, therefore it does not run Alpha
Wrap. The optional `--stage-dir` and `--component-csv` retain inspection
artifacts. `--preset local-diagnostic` is research-only (KNN=10, PCA K=20,
MAD=6) and is not the production default.

## Important source files

- `src/main.cpp`: command-line parsing, pipeline orchestration, CSV/PLY output,
  and `--filter-only` path.
- `src/denoise.cpp`, `include/denoise.hpp`: SOR, ROR, component clustering,
  robust-MAD local filtering, end metrics, and component metadata.
- `src/voxel_downsample.cpp`: centroid voxel downsampling.
- `src/alpha_wrap.cpp`: `CGAL::alpha_wrap_3` wrapper.
- `src/mesh_validation.cpp`: watertight/manifold/orientation/self-intersection
  and two-sided wrap detection.
- `src/volume.cpp`: signed tetrahedral volume computation in mm³.
- `filtering.py`: user-facing filtering-only CLI with research presets.
- `scripts/run_advanced_cleaning_benchmark.py`: process-isolated full experiment
  runner; `scripts/analyze_advanced_results.py` and
  `scripts/plot_advanced_cleaning.py` create summaries and figures.

## Formats and generated data

Input is binary little-endian PLY with XYZ; extra vertex properties such as RGB
are skipped by the C++ reader. Coordinates are millimetres. Output point clouds
are binary little-endian XYZ PLY. Mesh volume is mm³ / 1000 = mL.

Point clouds, meshes, benchmark CSVs, figures, CMake build products, and local
datasets are deliberately ignored by Git. Regenerate them locally.

## Current experimental findings

- Input used for the final residual analysis: 59,781-point SOR+ROR cloud.
- At 4 mm connectivity: main component 59,759 points; residual clusters 10, 6,
  5, and 1 points. Component filtering removes all 22 visible residual points.
- At 1.5 mm voxel on the dense cloud: 60,129 → 60,018 points; 111 off-body
  voxel representatives are removed.
- Main-component dominance is the safety condition: largest component must be
  at least 98% and second-largest below 1% of all points.
- Best volume candidate: voxel 1.5 mm, component radius 4 mm, alpha 6.25 mm,
  offset 0.10 mm; 388.686 mL; 1,359 vertices / 2,714 faces; watertight,
  manifold, one component, no two-sided wrap.
- 20-run best pipeline deployment time: mean 223.02 ms, P95 274.27 ms.
- Cluster-noise stress: component filtering estimated recovery 99.80%; local
  filtering alone was not sufficiently reliable.

## Assumptions and known limitations

- The recommended preset assumes a single cucumber is the dominant point-cloud
  component. For multiple objects or genuinely detached cucumber parts, inspect
  the component CSV and consider `--preset conservative` or
  `--preserve-detached-ends`.
- A 2.476 mm PCA-extreme change in the validated sample is entirely caused by
  visually confirmed detached noise; it is not evidence of cucumber shortening.
- There is no physical ground-truth volume, so reported volume differences and
  CV describe stability and geometric plausibility, not measurement accuracy.
- Alpha/offset choices remain data-dependent. Large alpha values can be
  topologically valid but geometrically over-smooth the cucumber.

## Build and validation

The tested environment is Windows + MSYS2 UCRT64 with CGAL, Eigen, TBB, OpenMP,
Ninja, GCC, Python, NumPy, and Matplotlib. Build Release with:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8
```

Run `python -m py_compile filtering.py` and a small `filtering.py` smoke test
after changing the CLI or filter-only path. There are currently no registered
CTest tests.

## Recommended next work

1. Validate the 1.5 mm / 4 mm component preset on independent cucumbers and
   cluttered captures.
2. Add labelled synthetic attached-spike data before reconsidering Adaptive
   Local filtering.
3. Obtain physical volume ground truth to measure absolute error.
4. Port and benchmark on Ubuntu using CGAL/TBB/OpenMP; preserve the two-sided
   wrap validation guard.
