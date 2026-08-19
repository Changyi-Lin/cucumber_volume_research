# Ring reconstruction + cubic B-spline endpoint benchmark

## Method and benchmark protocol

Release build; `std::chrono::steady_clock`; 10 warm-up runs and 100 measured runs per configuration. PLY read, PLY write, console output, and debug visualization are excluded. The same in-memory point cloud (1.5 mm voxel + 4.0 mm dominant component) is supplied to Alpha Wrap, Ring-flat, and Ring-spline. One-time preprocessing was 64.136 ms in this run.

Ring radii use exactly Q10/Q90 and `R=(Q10+Q90)/2`, with 48 angular bins and no radius calibration. `mean_ms`/percentiles include reconstruction + mesh validation + volume; `reconstruction_plus_volume_ms` excludes validation.

Reproduce with:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8 --target cucumber_ring_benchmark
.\build\cucumber_ring_benchmark.exe
```

## Question 1 — Ring count under 0.5 s and convergence

**Maximum tested valid Ring count under 0.5 s (p95): 148 Rings** (target 150, p95 39.008 ms).

Time-only capacity reached 296 actual Rings (target 300, p95 30.411 ms), but it is **not admissible** because the mesh self-intersects.

**Recommended Ring count: 90 Rings** (target 90). Its volume difference to the target-150 reference is 0.427%, while using the smallest tested topology-valid configuration whose reference error is at most 0.5% and remains converged at all higher topology-valid counts through the reference.

| target | actual | step mm | volume mL | change from previous % | diff to ref % | mean ms | p95 ms | watertight | manifold | self-intersection |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|
|20|20|13.508|264.334|0.000|7.405|13.147|14.116|yes|yes|no|
|30|30|9.005|274.278|3.762|3.922|14.899|16.014|yes|yes|no|
|40|40|6.754|278.090|1.390|2.586|16.593|17.244|yes|yes|no|
|50|50|5.403|280.303|0.796|1.811|17.644|18.701|yes|yes|no|
|60|60|4.503|281.661|0.484|1.335|19.700|20.655|yes|yes|no|
|75|75|3.602|283.331|0.593|0.750|22.544|23.493|yes|yes|no|
|90|90|3.002|284.255|0.326|0.427|24.326|25.406|yes|yes|no|
|120|120|2.251|285.700|0.508|0.079|31.734|33.875|yes|yes|no|
|150|148|1.801|285.473|0.079|0.000|37.236|39.008|yes|yes|no|
|180|178|1.501|286.457|0.345|0.345|21.614|22.969|yes|yes|yes|
|240|238|1.126|287.784|0.463|0.810|24.261|25.383|yes|yes|yes|
|300|296|0.901|287.896|0.039|0.849|28.400|30.411|yes|yes|yes|

Targets above the reference are retained as performance data but excluded from the geometric reference/recommendation because validation detected self-intersection or another topology failure. This is the observed high-resolution over-slicing limit.

Convergence thresholds (first tested count at or below threshold against reference):

- <=1.000%: 75 actual Rings
- <=0.500%: 90 actual Rings
- <=0.250%: 120 actual Rings
- <=0.100%: 120 actual Rings

## Question 2 — Endpoint B-spline cost

Recommended-count comparison (K=8, 4 generated rings per end):

| endpoint | fit ms | evaluation ms | ring generation ms | endpoint mesh ms | endpoint mean ms | endpoint p95 ms |
|---|---:|---:|---:|---:|---:|---:|
| Flat cap | 0 | 0 | 0.001|0.014|0.015| n/a |
| Cubic B-spline | 0.019|0.001|0.013|0.115|0.149|0.284|

The spline endpoint accounts for **1.479%** of Ring-spline reconstruction + volume time.

All endpoint sweeps are retained in `ring_benchmark.csv` (generated rings 2,3,4,5,6,8,10 at K=8; K=4,6,8,10,12 at four generated rings).

## Question 3 — Reconstruction comparison

Alpha Wrap parameters: alpha=6.250 mm, offset=0.100 mm.

| method | mean ms | p95 ms | volume mL | vertices | faces | watertight | manifold |
|---|---:|---:|---:|---:|---:|:---:|:---:|
|Alpha Wrap|122.698|130.408|388.686|1359|2714|yes|yes|
|Ring + flat cap|24.326|25.406|284.255|4322|8640|yes|yes|
|Ring + cubic B-spline tip|24.966|26.159|284.413|4706|9408|yes|yes|

All three methods meet 500 ms p95 on this machine. Ring-spline is the fastest topology-valid completed surface, but its volume is 26.827% below Alpha Wrap. Because there is no physical ground-truth volume, this benchmark establishes runtime and internal Ring convergence, not absolute accuracy or equivalence between reconstruction families. Alpha Wrap remains appropriate when continuity with the existing Alpha-derived volume is required.

Mean online core (one preprocessing pass + reconstruction + volume, validation excluded): Alpha Wrap 181.091 ms; Ring-flat 74.161 ms; Ring-spline 74.236 ms.

## Recommendation

Use `ring_spline` with target/actual 90/90 Rings, 48 angular bins, Q10/Q90 center radius, K=8 endpoint fit rings, and 4 generated rings per end. This is watertight/topology-valid, comfortably below 500 ms at p95, and uses the first volume-converged Ring count instead of spending the remaining budget on unnecessary axial resolution. The inspection mesh is `final_ring_spline_mesh.ply`.
