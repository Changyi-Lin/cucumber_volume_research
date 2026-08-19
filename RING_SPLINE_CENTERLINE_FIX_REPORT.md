# Ring spline reliable-boundary and centerline fix

## Root cause and geometry change

The radial template fix smoothed generated Rings, but the old robust mesh still retained unreliable outer Rings and extended only to the PCA point-cloud extreme. Their displaced centers created the boundary kink, while the approximately 0.24 radius-length extension produced a flat cap.

The new mesh selects the first reliable boundary from each endpoint, removes every Ring outside those boundaries from the side surface, fits a local endpoint centerline from five reliable body centers, and extrapolates along that fitted tangent/curvature for `tip_length_ratio * boundary_radius`. Radius uses a separate monotonic cubic Hermite taper; the previous four-Ring robust shape template remains unchanged.

Bottom reliable boundary: **Ring 2**; excluded Rings: **0–1**. Top reliable boundary: **Ring 86**; excluded Rings: **87–89**.

## Endpoint geometry metrics

| side | old boundary/turn max deg | new boundary tangent deg | new max turn deg | old L/R | new L/R | new tip length mm | monotonic radius |
|---|---:|---:|---:|---:|---:|---:|:---:|
| Bottom | 94.131979|0.861613|1.721669|0.246591|0.799970|10.067170|yes|
| Top | 39.094502|0.557367|1.114312|0.240542|0.799988|11.232343|yes|

The fitted path is position- and tangent-continuous at the reliable boundary. Generated local frames are recomputed from the fitted tangent, rather than inherited from a noisy outer Ring.

## Release benchmark

10 warm-up runs and 100 measured runs per configuration; `steady_clock`; I/O excluded. Total includes reconstruction, validation, and volume.

| method | mean ms | median ms | p95 ms | reconstruction ms | volume mL | vertices | faces | watertight | manifold | components | self-intersection | chi |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|---:|:---:|---:|
| Old robust endpoint | 43.155877|41.727600|53.756830|20.096726|285.719634|4706|9408|yes|yes|1|no|2|
| Centerline-fixed (ratio 0.8) | 40.925522|38.741300|57.335865|19.145862|285.411789|4466|8928|yes|yes|1|no|2|

Runtime change: **-5.168137%**. Volume change: **-0.107744%**. New p95 is **57.335865 ms**, well below 500 ms.

New mean endpoint costs: boundary selection 0.027525 ms; centerline fit 0.084334 ms; centerline evaluation 0.001699 ms; radius fit 0.000235 ms; Ring generation 0.038490 ms.

## Tip-length ratio sweep

| ratio | mean ms | median ms | p95 ms | volume mL | bottom max deg | top max deg | watertight |
|---:|---:|---:|---:|---:|---:|---:|:---:|
|0.500000|42.469721|39.230750|60.871880|283.865122|1.079134|0.697286|yes|
|0.750000|42.458950|39.582000|63.457610|285.156856|1.614980|1.044917|yes|
|0.800000|40.925522|38.741300|57.335865|285.411789|1.721669|1.114312|yes|
|1.000000|42.213624|39.219800|60.565250|286.420926|2.007618|1.391348|yes|

No ratio is calibrated to ground truth; 0.8 remains the default geometric prior.

## Reproduce

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8 --target cucumber_ring_robust_benchmark
.\build\cucumber_ring_robust_benchmark.exe
```
