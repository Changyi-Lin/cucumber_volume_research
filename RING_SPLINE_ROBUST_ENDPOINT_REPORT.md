# Robust Ring spline endpoint reconstruction

## Root cause

The old implementation copied all 48 radii from the outermost endpoint Ring and multiplied them by one spline scale. A common scalar cannot change coefficient of variation or normalized roughness, so any endpoint protrusion was propagated unchanged to every generated Ring. Increasing spline degree would only change the axial scale, not this shape error.

## Implementation

Each sector now applies median/MAD rejection (`3.0 * 1.4826 * MAD`) before Q10/Q90; MAD near zero, fewer than four samples, or fewer than three retained samples falls back safely. The 48 Q-center radii then use a circular five-sector median filter. Ring quality records point count, angular coverage, radius CV, and normalized circular roughness.

In each ten-Ring endpoint candidate region, shape eligibility requires robust relative CV/roughness thresholds (`median + 2.5 robust sigma`) capped at absolute 0.12 limits, plus robust lower bounds for coverage and point count. Radius fitting uses a deliberately looser coverage/point validity rule. Only the endpoint fitting sequence receives isotonic non-decreasing correction from tip toward body.

Bottom shape template Rings: **2, 3, 4, 5**. Top shape template Rings: **86, 85, 84, 83**. Each source Ring is normalized by its mean radius; a sector-wise median forms the template, which is normalized again to mean 1. The cubic B-spline supplies mean radius only. Tip shape regularization remains available but is disabled for this result.

The bottom template normalized-radius range is 0.951091 to 1.039394; the top range is 0.948633 to 1.037656. These non-zero ranges retain measured asymmetry rather than forcing a perfect circle. Endpoint rebound diagnostics flagged 0 candidate Rings in this dataset.

## Endpoint shape diagnostics

Direct re-analysis of the existing old mesh reproduces the reported endpoint pattern:

| old mesh Ring | mean radius mm | min mm | max mm | CV | normalized roughness |
|---:|---:|---:|---:|---:|---:|
|0|8.913139|6.200348|13.353423|0.184834|0.095224|
|1|10.073336|8.352365|11.282280|0.065618|0.064075|
|2|12.530414|11.052209|14.183089|0.053351|0.060718|
|87|11.921693|10.265959|13.634456|0.075171|0.040272|
|88|9.283848|7.293118|14.492353|0.125788|0.068994|
|89|6.378794|3.518881|10.169871|0.161277|0.146512|

| side/metric | old anchor | cleaned anchor | old generated mean | robust generated mean |
|---|---:|---:|---:|---:|
| Bottom CV | 0.184834|0.163915|0.184834|0.017021|
| Top CV | 0.161277|0.077331|0.161277|0.020258|
| Bottom normalized roughness | 0.095224|0.042681|0.095224|0.008717|
| Top normalized roughness | 0.146512|0.024239|0.146512|0.007516|

The robust generated Rings use a multi-Ring asymmetric template rather than a perfect circle; their CV/roughness no longer equals the noisy outer anchor. Inspect the template CSVs and original/cleaned endpoint PLYs for all 48 sectors.

## Release benchmark

10 warm-up runs and 100 measured runs per method; `std::chrono::steady_clock`; PLY I/O and console output excluded. Total includes reconstruction + validation + volume.

| method | mean ms | median ms | p95 ms | reconstruction ms | volume mL | watertight | manifold | components | self-intersection | chi | genus |
|---|---:|---:|---:|---:|---:|:---:|:---:|---:|:---:|---:|---:|
| Old spline | 40.479030|39.913750|45.101250|16.955084|284.413007|yes|yes|1|no|2|0.000000|
| Robust spline | 43.155877|41.727600|53.756830|20.096726|285.719634|yes|yes|1|no|2|0.000000|

Robust runtime change: **6.612923%**. Volume change: **0.459412%**. Robust p95 is 53.756830 ms, therefore it remains comfortably below 500 ms.

Robust-only mean costs: sector MAD 3.493351 ms; circular smoothing 0.561522 ms; Ring quality 0.021335 ms; shape template 0.027201 ms; endpoint radius fit 0.022569 ms; endpoint generation 0.026618 ms.

## Reproduce

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8 --target cucumber_ring_robust_benchmark
.\build\cucumber_ring_robust_benchmark.exe
```
