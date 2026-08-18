# 小黃瓜點雲 CGAL Alpha Wrap 3 實測結果

## Dataset

- 檔案：`pointcloud_3100.ply`（binary little-endian，Open3D 產生）
- 原始點數：1,573,158
- AABB：95.869 × 57.749 × 266.725 mm；對角線 289.253 mm
- PCA bounding box：53.532 × 87.989 × 272.794 mm
- 自訂 PLY reader 單次載入：約 0.52–0.59 s（Balanced 重複測試前為 549.3 ms）

PLY 載入屬 I/O／初始化成本，未混入題目定義的 online 核心時間。若每一 frame 都必須從硬碟重新讀完整 42.5 MB PLY，單是 I/O 就會超過 0.5 秒；實際即時系統應從相機／記憶體直接提供 XYZ buffer。

## Environment

- OS：Windows 11 Home 10.0.26200；本機沒有已安裝的 WSL，因此採 MSYS2 UCRT64
- CPU：Intel Core Ultra 7 256V，8 cores / 8 logical processors
- RAM：15.52 GiB
- CGAL：6.1-1，真正呼叫 `CGAL::alpha_wrap_3`
- Compiler：GCC 16.1.0；C++17
- CMake：4.3.3；TBB：2023.0.0
- Build：Release，`-O3 -march=native -DNDEBUG`

## Method 與有效性判定

流程為 binary PLY → packed-key centroid voxel → CGAL Alpha Wrap 3 → tiny-component cleanup → mesh validation → tetrahedral signed volume。

Alpha Wrap 在小 alpha 時會對少數孤立點產生 3–4 mm 的微小封閉 debris。僅在最大 raw component ≥95% faces 時保留最大 component；所有 raw component 數、移除數與 cleanup 時間都保留在 CSV。這不是用 convex hull、alpha shape 或 Poisson 取代 Alpha Wrap。

有效 Mesh 必須同時滿足：

- 單一 component、closed、triangle mesh、manifold、outward oriented
- 無 self-intersection
- Euler characteristic = 2、genus = 0（本報告三個推薦 mesh 均如此）
- two-sided 偵測未觸發

two-sided 偵測在 PCA 主軸 20%、35%、50%、65%、80% 位置，通過局部點雲中心沿兩個橫向各投一條無限直線，共 10 條；交點去重後，若中位數 ≥4 且至少 40% 線有 ≥4 個交點，即標記 `INVALID_TWO_SIDED_WRAP`。中央五個實際截面圖亦用精確 triangle-plane intersection 驗證外殼只有單一輪廓。

## Voxel downsampling benchmark

| Voxel (mm) | 點數 | 單次時間 (ms) |
|---:|---:|---:|
| Original | 1,573,158 | 0 |
| 0.50 | 736,676 | 396.71 |
| 0.75 | 322,475 | 78.25 |
| 1.00 | 164,821 | 46.13 |
| 1.25 | 95,054 | 30.43 |
| 1.50 | 60,129 | 20.84 |
| 1.75 | 40,797 | 16.07 |
| 2.00 | 29,223 | 13.61 |
| 2.50 | 16,760 | 14.72 |
| 3.00 | 10,745 | 15.18 |

這是各 voxel 第一次、各一次的獨立 benchmark；Balanced 的 20 次 downsample 平均 41.63 ms、中位數 35.70 ms、P95 73.37 ms。

## Parameter search

最終 `alpha_wrap_benchmark.csv` 有 270 筆：

- 規定 coarse grid：voxel 1/1.5/2/2.5/3；alpha 20/25/30/35/40；offset 0.2/0.3/0.5/0.75
- 初始完整 alpha：15/20/25/30/35/40/45/50 mm
- offset：0.10/0.20/0.30/0.50/0.75/1.00 mm
- fine grid：voxel 1.25–2.5、alpha 5/6.25/7.5/8.75/10/12.5/15/17.5、offset 0.2/0.3/0.5
- 另補 voxel 0/0.5/0.75/1.0 與 offset 0.1 的品質／即時候選，以及 alpha 3 的明確失敗案例

主要規律：

- alpha 變大會跨過彎曲凹部並明顯放大體積。例如 voxel 2、offset 0.3 下，alpha 15→50 的體積約 421.1→505.0 mL。
- offset 變大符合 conservative wrap 預期地放大體積。voxel 1.5、alpha 15 下，offset 0.1→1.0 的體積 416.82→447.05 mL，增加 30.23 mL（7.25%）。
- alpha 5 位於「封住開口」與「鑽入孔洞」的臨界區，附近參數不連續；alpha 6.25 起在測試區內沒有 two-sided case。

## Two-sided wrapping 與 failure cases

共找到 5 個 `INVALID_TWO_SIDED_WRAP`：

| Voxel | Alpha | Offset | 錯誤體積 (mL) | 判定 |
|---:|---:|---:|---:|---|
| 2.0 | 3.0 | 0.3 | 170.41 | 10/10 中央線為四交點 |
| 1.5 | 5.0 | 0.3 | 218.38 | 中位四交點 |
| 1.75 | 5.0 | 0.2 | 210.13 | 中位四交點 |
| 2.25 | 5.0 | 0.2 | 201.61 | 中位四交點 |
| 2.5 | 5.0 | 0.5 | 210.95 | 中位四交點 |

它們在 cleanup 後仍可 technically closed/manifold；若只檢查 `is_closed=true` 就會把約 170–218 mL 的薄表面殼誤當小黃瓜實體。alpha 5 的其他組合可得 375–396 mL 的正常單殼，正是這種參數不連續性使其不適合當量產 Balanced 設定。

## Volume stability

穩定區選定 alpha 6.25、offset 0.10，voxel 1.0/1.25/1.5/1.75/2.0：

```text
393.866, 391.429, 389.410, 386.994, 384.974 mL
```

- mean：389.335 mL
- standard deviation：3.144 mL
- CV：0.807%
- max-min：8.892 mL，即 mean 的 2.284%

這是一個隨 voxel 平滑變化的區域，而不是 200↔390 mL 的跳變。沒有實際排水體積 Ground Truth，因此不能宣稱其中哪一個數字「最準」；以下 Best Quality 與 Balanced 是貼合度、穩定性與速度的工程選擇。

## 三種推薦設定與重複 runtime

| 方案 | Voxel | 點數 | Alpha | Offset | Volume | Vertices / Faces | Online mean / median / P95 | 判定 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Fastest | 3.0 | 10,745 | 40 | 0.75 | 517.063 mL | 135 / 266 | 20.28 / 20.15 / 22.06 ms | valid，但幾何很粗且體積偏大 |
| Best Balanced | 1.0 | 164,821 | 6.25 | 0.10 | **393.866 mL** | 1,515 / 3,026 | **262.40 / 248.68 / 337.14 ms** | valid，最推薦 |
| Best Quality | Original | 1,573,158 | 5.0 | 0.10 | 391.253 mL | 2,348 / 4,692 | 1166.14 / 1161.66 / 1254.15 ms | valid，超過 0.5 s |

Fastest 與 Best Quality 各 10 次，Balanced 20 次；每組另做一次 warm-up 並排除。三組每次結果均為 watertight、manifold、無 self-intersection、單一 component、無 two-sided wrap。

### Balanced 20 次詳細時間（ms）

| 階段 | Mean | Median | Std | Min | Max | P95 |
|---|---:|---:|---:|---:|---:|---:|
| Downsample | 41.63 | 35.70 | 13.67 | 31.63 | 83.15 | 73.37 |
| Alpha Wrap | 220.27 | 212.00 | 19.47 | 203.44 | 277.93 | 265.47 |
| Cleanup | 0.484 | 0.450 | 0.125 | 0.316 | 0.841 | 0.667 |
| Volume | 0.016 | 0.015 | 0.004 | 0.014 | 0.029 | 0.023 |
| Reconstruction only (Alpha + Volume) | 220.28 | 212.02 | 19.47 | 203.46 | 277.95 | 265.48 |
| Online core | **262.40** | **248.68** | **31.62** | **239.04** | **338.16** | **337.14** |
| Validation | 13.31 | 12.97 | 0.94 | 11.60 | 15.74 | 14.58 |
| Total including validation | 275.71 | 261.28 | 32.29 | 252.01 | 352.83 | 352.69 |

Balanced 的 Alpha Wrap 約占 online mean 84%，是主要 bottleneck；downsample 約 16%；volume 本身只有約 0.016 ms。

## Number of points vs runtime

固定 alpha 15、offset 0.3，從原始點做 deterministic systematic sample；原始點列不做 copy/downsample：

| Input points | Alpha Wrap (ms) | Online core (ms) |
|---:|---:|---:|
| 10,000 | 27.95 | 28.27 |
| 20,000 | 36.49 | 36.79 |
| 30,000 | 41.48 | 42.20 |
| 50,000 | 55.95 | 56.69 |
| 75,000 | 66.97 | 67.84 |
| 100,000 | 78.40 | 79.90 |
| 150,000 | 99.47 | 102.16 |
| 300,000 | 166.38 | 169.64 |
| 500,000 | 244.80 | 248.5 |
| 750,000 | 360.1 | 364.8 |
| 1,000,000 | 489.9 | 496.3 |
| 1,250,000 | 581.9 | 588.5 |
| 1,573,158 | **742.0** | **742.0** |

單次量測的 500 ms 交界落在 1.0M–1.25M 點；1.0M 的 496.3 ms 幾乎沒有抖動餘裕，工程上要「穩定」低於 0.5 s 建議 ≤750k 點。推薦 Balanced 直接降到 164,821 點，其 20 次 P95 337.1 ms，安全餘裕大得多。

原始 157 萬點的時間取決於幾何參數：alpha 15 / offset 0.3 為 742.0 ms；較精細的 alpha 6.25 / offset 0.1 單次 895.9 ms；Best Quality alpha 5 / offset 0.1 的 10 次 Alpha Wrap mean 1165.4 ms、P95 1253.1 ms。

## Figures 與 Mesh

`results/figures/` 包含：

- `runtime_vs_voxel_size.png`
- `volume_vs_voxel_size.png`
- `volume_vs_alpha.png`
- `volume_vs_offset.png`
- `runtime_vs_volume.png`
- `runtime_volume_stability_pareto.png`
- `points_vs_alpha_wrap_runtime.png`
- 原始點雲、Balanced mesh、overlay、top/bottom local end-plane、五個 PCA cross-sections

`results/meshes/` 內的 `best_fast_mesh`、`best_balanced_mesh`、`best_quality_mesh` 均同時輸出 binary PLY 與 STL。

## 對 12 個問題的直接回答

1. **157 萬原始點直接做 Alpha Wrap：** alpha 15 / offset 0.3 實測 741.95 ms；品質設定會增至約 0.90–1.17 s。
2. **多少點可進 0.5 s：** 單次臨界約 1.0M；要有穩定餘裕建議 ≤750k。Balanced 的 164,821 點 P95 337.1 ms。
3. **最合適 voxel：** 1.0 mm；保留 164,821 點，P95 遠低於 500 ms，且比 1.5–2.0 mm 保留更多細節。
4. **最穩定封兩端的 alpha：** 6.25 mm；5 mm 在部分鄰近參數會鑽洞，15 mm 以上雖穩定但體積膨脹較明顯。
5. **offset：** 0.10 mm；在 alpha 6.25 區域穩定、複雜度可控，且是測試中 conservative expansion 最小者。
6. **two-sided：** 有，集中在 alpha 3 與 alpha 5 臨界區，完整組合列於上表。
7. **最佳 Alpha Wrap 體積：** 推薦 Balanced 為 **393.866 mL**；Best Quality 原始點設定為 391.253 mL。無 Ground Truth，不宣稱前者是真實體積。
8. **合理參數體積變異：** alpha 6.25 / offset 0.1 / voxel 1–2 mm 的 CV 0.807%，range 8.892 mL（2.284%）。
9. **Alpha Wrap + Volume <0.5 s：** 是；Balanced reconstruction-only P95 265.48 ms，含 downsample/cleanup 的 online P95 337.14 ms。
10. **不能時的 bottleneck：** 原始點／細參數時 Alpha Wrap 是主因；PLY I/O 若每 frame 重讀亦約 0.55 s，但應視為可避免的初始化／資料介面成本。
11. **Best Balanced：** voxel 1.0 mm、164,821 點、alpha 6.25 mm、offset 0.10 mm。
12. **是否適合大量即時估計：** 對此資料，在記憶體輸入且先 voxel 到約 165k 點時適合 CPU 0.5 s 級處理；必須固定 two-sided/拓樸 guard。若目標是更高 frame rate（例如 10–30 FPS），Alpha Wrap 仍太慢，應評估更輕量的切片／高度圖或 GPU 方法。

## Final conclusion

```text
Recommended configuration (Best Balanced)

Voxel size:       1.0 mm
Input points:     164,821
Alpha:            6.25 mm
Offset:           0.10 mm

Downsample:       mean 41.63 ms, P95 73.37 ms
Alpha Wrap:       mean 220.27 ms, P95 265.47 ms
Cleanup:          mean 0.484 ms
Volume:           mean 0.016 ms
Online total:     mean 262.40 ms, median 248.68 ms, P95 337.14 ms

Volume:           393.866 mL
Watertight:       yes
Manifold:         yes
Self-intersection:no
Two-sided wrap:   no

Conclusion: CGAL Alpha Wrap + Volume can reliably meet 0.5 s on this dataset
after 1.0 mm voxel downsampling. Direct 1.57M-point processing cannot reliably
meet 0.5 s, and alpha <= 5 mm is not a stable production region because some
nearby settings wrap both sides of the open surface.
```

## Denoising Study（補充實驗）

### 實驗設計與執行情形

在既有 benchmark 上加入 `Voxel → SOR / ROR / SOR+ROR → Alpha Wrap`。主 CSV 共 **348 次實際執行**：65 組 coarse、53 組細搜、12 組端部保護、144 組 Alpha/offset/voxel 鄰域比較、6 組 clustering、64 組人工 outlier stress、4 組代表點雲輸出；另有四方法各 10 次、共 40 次交錯重複測量。

- SOR 涵蓋題目列出的 K = 10/20/30/50/80 與 std ratio = 0.75/1/1.25/1.5/2/2.5/3；先 coarse 再補邊界值，沒有做無意義的完整笛卡兒積。
- ROR 涵蓋 radius 1–5 mm、min neighbors 2/3/5/8/10/15；正式比較採 adaptive radius = `2.5 × voxel size`。
- 端部保護比較 0/5/7.5/10%；以 PCA 主軸投影定義兩端。
- Alpha Wrap 比較 voxel 1/1.5/2/2.5 mm、alpha 6.25/7.5/10 mm、offset 0.1/0.2/0.3 mm。前一輪已證明 alpha 20–40 mm 會過度跨越彎曲凹部，因此本輪依「前一輪合理區間」使用 6.25–10 mm，而不重複把已知不合理區域當成候選。

所有 SOR/ROR 查詢使用 CGAL Kd-tree；逐點鄰域查詢以 OpenMP 8 threads 執行。`deployment_total_ms` 包含 voxel、denoise、Alpha Wrap、必要 component cleanup、volume；`research_total_ms` 再加入 mesh validation。

### 去噪參數與端部保護選擇

未保護的 12 組 coarse SOR 全部觸發 `POSSIBLE_OVER_DENOISE`。選定的 SOR K=30/std=3.0 在 global 模式會讓 PCA 長度縮短 **3.382 mm（1.241%）**，底部刪除 1.677%；保護兩端各 5% 後，長度變化與頂／底刪除率都成為 0。

| Method | End protection | Removed | Top / middle / bottom removed | Length loss | Valid |
|---|---:|---:|---:|---:|---|
| SOR K30/std3 | 0% | 0.640% | 0.627 / 0.611 / 1.677% | 3.382 mm | No：過度去噪 |
| SOR K30/std3 | 5% | 0.579% | 0 / 0.611 / 0% | 0 mm | Yes |
| ROR r3.75/n3 | 0% | 0.0067% | 0.063 / 0.0035 / 0.062% | 0 mm | Yes |
| ROR r3.75/n3 | 5% | 0.0033% | 0 / 0.0035 / 0% | 0 mm | Yes |
| SOR+ROR | 5% | 0.579% | 0 / 0.611 / 0% | 0 mm | Yes |

5%、7.5%、10% 都保住長度；5% 對中段保留最多可用 filtering 範圍，因此採 5%。沿兩端建立 0–2、2–4、4–6、6–8、8–10 mm 截面後，四方法的 robust width、height、ellipse area 相對 voxel-only 都是 **0% 變化**，且這 0–10 mm 內沒有刪除點。

### 自然點雲的局部密度與 outlier 情形

voxel 1.5 mm 將 1,573,158 點聚合為 60,129 點，數量下降 96.18%；這是空間聚合，不應等同「刪除了 96.18% noise」。在 voxel-only 點雲上，溫和 ROR 只找到 **2 個**孤立點（0.0033%），表示本資料降採樣後沒有大量明顯離群點；SOR 會額外刪 348 個（0.579%）中段統計尾端點。

Voxel-only 最近鄰距離中位數為：bottom 1.097 mm、middle 1.109 mm、top 1.073 mm。就局部 NN spacing 而言，頭尾並沒有比中段更稀疏；頭尾總點數較少主要也與截面面積縮小有關。不過 global SOR 的實測長度損失仍證明端部保護有必要，不能只靠全域密度假設。

### 四方法最終比較

下表固定 voxel 1.5 mm、alpha 6.25 mm、offset 0.1 mm；runtime 是 10 次交錯重複平均。Volume CV 是每個 voxel 內對 3 alpha × 3 offset 計算 CV，再對四個 voxel 取平均；沒有 Ground Truth，因此只能比較穩定性，不能宣稱 accuracy 改善。

| Metric | No denoise | SOR | ROR | SOR+ROR |
|---|---:|---:|---:|---:|
| Remaining points | 60,129 | 59,781 | 60,127 | 59,781 |
| Removed % | 0 | 0.579 | 0.0033 | 0.579 |
| End damage / length loss | 0 mm | 0 mm | 0 mm | 0 mm |
| Denoise mean (ms) | 0 | 131.73 | 44.40 | 165.01 |
| Alpha Wrap mean (ms) | 142.58 | 141.64 | 138.56 | 131.29 |
| Deployment mean / P95 (ms) | 160.88 / 169.89 | 291.20 / 329.99 | **201.15 / 216.82** | 314.27 / 323.44 |
| Volume (mL) | 389.410 | 387.648 | **389.346** | 387.648 |
| Mean within-voxel volume CV | 1.843% | 1.599% | 1.814% | **1.597%** |
| Vertices / faces | 1,440 / 2,876 | 1,399 / 2,794 | 1,423 / 2,842 | 1,399 / 2,794 |
| Watertight / manifold | Yes / Yes | Yes / Yes | Yes / Yes | Yes / Yes |
| Two-sided wrap | No | No | No | No |
| Recommendation | 最快；受控乾淨輸入 | 穩定性較佳但成本高 | **Best overall** | CV 最低但重複過濾不划算 |

在自然資料上，ROR 只讓 Alpha Wrap mean 從 142.58 降至 138.56 ms（約 2.8%），但完整 total 增加約 40.27 ms；因此不能說自然資料靠去噪加速了整條 pipeline。SOR/組合把 faces 降約 2.9%，CV 相對下降約 13%，但 volume 也改變約 −1.76 mL；沒有 GT，這不是 accuracy 證據。

### 孤立群集檢查

radius 3.0/3.75/4.5 mm 找到 12–14 群，最大群約 60,011–60,018 點、第二群只有 17 點。若不保護端點而直接留最大群，會移除 111–118 點並讓長度縮短 **2.476 mm**，被判定過度去噪。保護兩端各 5% 後只移除 89–93 點，長度損失 0。由於自然 ROR 已只刪 2 點、最終 mesh 也正常，預設 pipeline 不啟用 clustering；它保留為診斷／明顯小群污染時的選項。

### 人工 outlier stress test

在 voxel 1.5 mm 點雲加入 0.01/0.05/0.1/0.5/1% 離群點，置於 AABB 外 5/10/20 mm；四方法共 64 次（含 clean baseline）。

| Method | Valid stress cases | Mean abs volume change | Worst change | 1% / 20 mm 結果 |
|---|---:|---:|---:|---|
| No denoise | 4/15 (26.7%) | 1.898% | 14.190% | 0 點刪除；393.973 mL；multiple components + two-sided |
| SOR | 15/15 | 0.421% | 0.723% | 636 點刪除；389.580 mL；valid |
| ROR | **15/15** | **0.000%** | **0.000%** | 603 點刪除；389.346 mL；valid；198.43 ms |
| SOR+ROR | 15/15 | 0.339% | 0.413% | 638 點刪除；389.174 mL；valid |

這是 ROR 最有價值的結果：乾淨點雲幾乎不動，加入離群點後則把污染刪掉，15 組 volume 全部回到各自 clean baseline，且保持單一 watertight/manifold、無 two-sided wrap。Voxel-only 在 1%/5 mm 時體積膨脹 14.19%，也常因額外 component 或 two-sided 判定失敗。

### 對 Denoising Study 11 個問題的回答

1. **是否有明顯自然 outlier：** voxel 1.5 mm 後只有 2 個點被保守 ROR 判為孤立，沒有大量自然 outlier。
2. **Voxel 本身移除多少 noise：** 點數聚合下降 96.18%，但 voxel 不是 outlier classifier，不能把這個比例稱作 noise；兩個自然孤立點仍會保留。
3. **SOR 移除多少：** 選定 protected SOR 移除 348 點、0.579%，全在 middle metric region。
4. **ROR 移除多少：** 自然資料 2 點、0.0033%；1%/20 mm stress 時移除 603 點。
5. **被刪除點在哪裡：** SOR 主要是中段／保護區外表面的統計尾端；ROR 是兩個孤立點。PLY 與全體／端部放大圖均已輸出。
6. **頭尾是否誤刪：** global SOR 會；5% protection 後 0–10 mm 截面、端點點數與長度均無損。
7. **端部最安全方法：** protected ROR；自然資料只刪 2 點且端點 0 點。
8. **是否改善 Alpha Wrap runtime：** ROR 約快 2.8%，SOR+ROR 約快 7.9%；但 denoise 成本使完整 total 都比 voxel-only 慢。
9. **是否改善 volume stability：** SOR/組合的 CV 明顯較低；ROR 僅小幅改善。無 GT，不能推論 accuracy。
10. **是否仍 <0.5 s：** voxel 1.5 mm 的四方法 mean/P95 全低於 500 ms；ROR P95 216.82 ms。voxel 1.0 mm 的單次 SOR/組合曾超過 500 ms，因此不作去噪版最佳部署設定。
11. **是否值得加入最終 pipeline：** 若輸入受控且保證無 stray points，voxel-only 最快；若實際農場／相機流程可能出現孤立點，protected adaptive ROR 值得加入，stress robustness 的提升遠大於約 44 ms 成本。

### 分類推薦

- **Best Conservative Denoising：** ROR，radius = 2.5 × voxel、min neighbors = 3、end protection = 5%。
- **Best Runtime Denoising：** ROR；若允許完全不去噪，voxel-only 仍最快。
- **Best Stable Volume Denoising：** SOR+ROR 的 1.597% 最低，但僅比 SOR 1.599% 好 0.002 percentage point；工程上選 SOR 較合理。
- **Best Overall：** protected adaptive ROR；自然資料低侵入、stress 15/15 有效、P95 216.82 ms。

### 補充實驗最終設定

```text
Recommended preprocessing

Voxel size:              1.5 mm
Denoising:               Radius Outlier Removal (adaptive)
SOR K / std ratio:       N/A
ROR radius:              3.75 mm = 2.5 × voxel
ROR min neighbors:       3
End protection:          top 5% + bottom 5% of PCA length
Clustering:              disabled by default

Original points:         1,573,158
After voxel:             60,129
Remaining after ROR:     60,127
Denoise removed:         2 points (0.0033%)

Voxel runtime:           mean 17.79 ms
Denoise runtime:         mean 44.40 ms
Alpha Wrap runtime:      mean 138.56 ms
Cleanup runtime:         mean 0.39 ms
Volume runtime:          mean 0.014 ms
Deployment total:        mean 201.15 ms, P95 216.82 ms
Research total:          mean 209.86 ms

Volume:                  389.346 mL
Volume stability:        mean within-voxel CV 1.814%

Watertight:              Yes
Manifold:                Yes
Two-sided wrap:          No
End geometry preserved:  Yes

Does denoising improve Alpha Wrap?
Natural clean data:      Marginal (slightly lower wrap time/CV, higher total time)
Contaminated data:       Yes (15/15 valid; 0% stress volume drift for ROR)

Recommended for real-time deployment? Yes
```

所有原始列、衍生表、截面量測、密度統計、剩餘／刪除 PLY 與圖表分別位於 `results/denoising_benchmark.csv`、`denoising_repeatability.csv`、`denoising_method_summary.csv`、`outlier_stress_*.csv`、`end_cross_section_metrics.csv`、`density_analysis.csv`、`denoised_points/` 與 `figures/`。

## Residual Noise Analysis（進階群集去噪，取代上一節的 clustering 結論）

### 資料與判定範圍

本輪實際使用專案內 `results/denoised_points/sor_ror_remaining.ply`；它恰為 **59,781 點**，且在 4 mm 連通條件下重現題目所列的 `59,759 + 10 + 6 + 5 + 1` 結構。題名含 `(1)` 的檔案並未出現在 workspace，因此不可宣稱另讀到一個未提供的不同檔案。

| Dataset | Points | AABB (mm) | PCA extents (mm，小→大) |
|---|---:|---:|---:|
| Raw `pointcloud_3100.ply` | 1,573,158 | 95.869 × 57.749 × 266.724 | 53.532 × 87.989 × 272.794 |
| SOR + ROR remaining | 59,781 | 83.186 × 56.938 × 265.965 | 50.760 × 64.163 × 272.626 |

10% PCA 分區的 nearest-neighbor 統計如下。原始密集點雲 overall mean/median 為 0.275/0.292 mm；SOR+ROR 後為 1.090/1.109 mm。後者 bottom/middle/top median 分別為 1.099/1.111/1.101 mm，沒有足夠差異支持「端部稀疏點一律是 noise」。5%、7.5%、10% 的完整 mean、median、STD、P5/P25/P75/P95 在 `advanced_density_regional_summary.csv`。

### 本輪實作

Release C++ pipeline 現在支援：

```text
Voxel
→ SOR / ROR（可選）
→ radius connected components
→ regional robust-MAD local density
→ local PCA plane residual + surface variation
→ weighted outlier score
→ real CGAL::alpha_wrap_3
→ topology/two-sided validation
→ volume
```

Component metadata 包含 ID、點數、centroid、AABB、至主體最近距離、是否保留；不同群集輸出 RGB PLY。Adaptive Local 使用 density/plane/surface variation 權重 0.5/0.3/0.2，且必須同時滿足低密度與幾何不一致條件才刪除。程式另記錄已人工確認之 disconnected-noise length contribution，保留 raw length change，同時把真正無法解釋的幾何縮短獨立為 `unexplained_length_change_mm`。

### 為什麼 SOR + ROR 後仍看得到 noise

SOR 判斷的是相對全域／區域統計尾端，ROR 判斷的是半徑內鄰居數。10、6、5 點的小群內部彼此有足夠鄰居，因此不是典型 single-point sparse outlier；它們雖與小黃瓜主體分離，仍可同時通過 SOR 與 ROR。這是 **cluster-type noise**，不是把 SOR/ROR threshold 再調強就能安全解決的問題。

### Connected Component 結果

在 residual cloud 上，3.5–5.0 mm 都得到相同的五群；採 **4.0 mm**，因為它位於穩定 plateau 中央，亦等於 1.5 mm voxel 的 2.67 倍。

| Component | Points | Distance to main (mm) | Visual result |
|---|---:|---:|---|
| Main cucumber | 59,759 | 0 | 保留 |
| Noise 1 | 10 | 8.788 | 刪除 |
| Noise 2 | 6 | 5.358 | 刪除 |
| Noise 3 | 5 | 9.308 | 刪除 |
| Noise 4 | 1 | 5.356 | 刪除 |

最大群占 99.963%，第二大群僅 0.0167%；22 點全部在共同 PCA 視角中位於主體之外。原始點雲經 1.5 mm voxel 後則為主體 60,018 點，加上 11 個小群（17、17、13、13、13、10、9、7、6、5、1），共刪 111 點；這些點也全部在 `final_component_removed.png` 中位於表面之外。

raw PCA extreme length 的確縮短 2.476 mm，但這 2.476 mm **全由已確認的離體群集貢獻**；`unexplained_length_change_mm = 0`。這符合題目所允許的人工檢查例外，不能再把它誤報為小黃瓜端部被削掉。

### Adaptive Local 是否必要

保守掃描結果中，MAD 3–5 會刪 331–2,405 個點並造成明顯幾何縮短，已停止往更強方向測。最保守可用的診斷設定為 KNN K=10、local PCA K=20、regional MAD=6、top/bottom threshold ×3：在 component 後再標出 133 點，且沒有再刪 top/bottom 10 mm 的點。

但共同視角圖顯示這 133 點廣泛散布在正常外表面，沒有形成可辨識毛刺；cluster stress 中 Adaptive Local 單獨使用的平均估計恢復率只有 0.77%，且僅 9/20 case 通過完整 validation。Local PCA／plane residual／surface variation 並未提供足以抵銷誤刪風險的改善。因此：

- MAD=6 是「若未來新資料真的出現 attached spike」時的保守診斷起點。
- 本資料的 production pipeline **關閉 Adaptive Local**。
- 沒有證據顯示仍存在貼在表面的真 noise；不能把正常曲率／取樣變化硬稱為 noise。

### 題目指定的 residual 直接比較

| Stage | Remaining | Removed | Raw length change | Confirmed-noise length | Unexplained damage | Top / middle / bottom removed | Volume |
|---|---:|---:|---:|---:|---:|---:|---:|
| SOR + ROR before | 59,781 | 0 | 0 | 0 | 0 | 0 / 0 / 0% | 387.648 mL |
| After Component | 59,759 | 22 | 2.476 mm | 2.476 mm | **0 mm** | 0.022 / 0 / 0.453% | 386.865 mL |
| After Adaptive Local | 59,626 | 155 total | 2.476 mm | 2.476 mm | **0 mm** | 0.022 / 0.262 / 0.453% | 385.777 mL |

Component 使 residual volume 改變 −0.783 mL（−0.202%）；Local 再改變 −1.088 mL。無 Ground Truth，這些只能描述敏感度，不能稱為 accuracy 改善。

端部 0–2、2–4、4–6、6–8、8–10 mm 截面使用確認後的 main-component endpoints 對齊。Component 在 bottom 前三 bin 刪除 5/3/2 個已確認 noise，top 第一 bin 刪 1 個 noise，其餘 bin 為 0；Adaptive Local 在兩端 0–10 mm **沒有再刪任何點**。因此截面改變源自移除離體點，不是主體寬、高或面積被削除。

### 七組 Denoising Ablation（voxel 1.5 mm）

固定 alpha=6.25 mm、offset=0.10 mm。表中 End Damage 是扣除人工確認 component-noise contribution 後的 unexplained damage。

| Pipeline | Remaining Points | Removed % | End Damage | Cleaning Time | Alpha Wrap Time | Total | Volume | Valid |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Voxel | 60,129 | 0 | 0 mm | 0.00 ms | 129.94 ms | 146.89 ms | 389.410 mL | Yes |
| Voxel + SOR | 59,822 | 0.511 | 0 mm | 43.94 ms | 119.74 ms | 178.75 ms | 387.607 mL | Yes |
| Voxel + ROR | 60,127 | 0.0033 | 0 mm | 22.41 ms | 130.37 ms | 167.22 ms | 389.346 mL | Yes |
| Voxel + SOR + ROR | 59,822 | 0.511 | 0 mm | 57.39 ms | 120.51 ms | 192.10 ms | 387.607 mL | Yes |
| + Component | 59,800 | 0.547 | 0 mm | 98.07 ms | 117.85 ms | 232.23 ms | 387.294 mL | Yes |
| + Adaptive Local | 59,668 | 0.767 | 0 mm | 148.41 ms | 119.03 ms | 282.50 ms | 386.148 mL | Yes |
| Component + Adaptive Local | 59,890 | 0.397 | 0 mm | 103.12 ms | 123.08 ms | 242.25 ms | 387.221 mL | Yes |

額外的最終候選 **Voxel + Component（不含 Local）** 為 60,018 點。它清除實際 cluster noise，少了 SOR 的 307 個正常統計尾端點，也少了 Local 的 128 個正常表面候選。因此 SOR/ROR 並非此資料最合理的 production 前處理。

### Cluster-noise Stress Test

在 voxel cloud 加入 clustered noise 0.01/0.05/0.1/0.5/1%，距表面 3/5/10/20 mm，共四方法 × 20 個污染 case。

| Method | Mean estimated recovery | Mean abs volume delta | Worst delta | Valid cases |
|---|---:|---:|---:|---:|
| SOR + ROR | 27.12% | 1.057 mL | 1.493 mL | 20/20 |
| Component | **99.80%** | **0.012 mL** | **0.119 mL** | **20/20** |
| Adaptive Local | 0.77% | 1.719 mL | 9.151 mL | 9/20 |
| Component + Adaptive Local | 99.80% | 0.022 mL | 0.219 mL | 20/20 |

Component 在最難的 1%/3 mm case 恢復 593/601 點（98.67%），5/10/20 mm 都為 601/601；附加 Local 沒有提高 recovery，反而增加 baseline surface removal 與 volume drift。

### Cleaning 後 Alpha Wrap、Mesh 與 Volume

新的 component-cleaned cloud 完整測了 alpha 20/25/30/35/40/45 × offset 0.1/0.2/0.3/0.5/0.75/1.0，以及 fine alpha 6.25/7.5/10/15/17.5 × offset 0.1/0.2/0.3/0.5。56/56 均為單一、closed、watertight、manifold、outward、無 self-intersection、無 two-sided wrap。

不過 alpha 20–45 把體積推高到約 412–455 mL，顯示跨越彎曲凹部；有效拓樸不代表幾何合理。最保真候選仍為 alpha=6.25、offset=0.10 mm。

- voxel-only → component-cleaned：volume 389.410 → 388.686 mL（−0.725 mL，−0.186%）。
- vertices/faces：1,440/2,876 → 1,359/2,714（faces −5.63%）。
- 同一組 Alpha Wrap 單次約 129.94 → 119.64 ms；完整 pipeline 因 component search 多出成本，total 並未變快。
- 合理鄰域 alpha 6.25/7.5、offset 0.1/0.2/0.3 的 volume CV：voxel-only 0.986%，component-cleaned **0.867%**，屬小幅改善。
- 固定設定 20 次的 volume 完全相同（CV 約 0）；這只證明實作 deterministic，不代表測量 accuracy。

三個輸出 mesh：

| Candidate | Alpha / offset | Volume | Vertices / faces | Single-run total | Comment |
|---|---:|---:|---:|---:|---|
| Best Quality / Overall | 6.25 / 0.10 mm | **388.686 mL** | 1,359 / 2,714 | 182.52 ms | 最保真、正式推薦 |
| Best Balanced | 7.5 / 0.20 mm | 396.970 mL | 863 / 1,722 | 144.09 ms | 較平滑；體積增加 |
| Best Fast | 10 / 0.30 mm | 407.496 mL | 541 / 1,078 | 123.09 ms | 僅作速度候選，不是體積首選 |

### 20 次正式 Runtime（Best Overall）

第一輪 warm-up 排除；PLY load 不含在 deployment runtime。

| Metric | Mean | Median | STD | Min | Max | P95 |
|---|---:|---:|---:|---:|---:|---:|
| Voxel | 21.85 | 16.26 | 7.55 | 13.97 | 33.15 | 32.14 ms |
| Component | 56.90 | 47.22 | 14.13 | 40.62 | 75.09 | 74.19 ms |
| Alpha Wrap | 141.40 | 141.79 | 21.33 | 118.01 | 166.83 | 166.62 ms |
| Volume | 0.014 | 0.014 | 0.0012 | 0.0125 | 0.0179 | 0.0161 ms |
| Deployment total | **223.02** | **210.31** | 41.79 | 177.09 | 275.46 | **274.27 ms** |

Alpha Wrap 約占 mean total 63%，是最大 bottleneck；component 約 26%，voxel 約 10%。20/20 都明顯小於 500 ms。

### 本輪 20 個問題的直接回答

1. SOR+ROR 殘留是因為小群內點互相支持，不符合 single sparse outlier 假設。
2. Connected Component 有效；residual 的 22/22 個離體點全清除。
3. 最佳 radius 為 4.0 mm；3.5–5.0 mm 是相同輸出的穩定 plateau。
4. Adaptive Local 對本資料沒有必要。
5. 若作診斷，最佳保守 MAD factor 為 6；production 關閉。
6. 沒有找到可被可信區分的 attached surface noise。
7. Local surface consistency 沒有改善可見 noise，反而標出正常表面點。
8. residual raw length 差 2.476 mm，全部是已確認 noise；unexplained geometry loss 0 mm。
9. 最終 raw-voxel pipeline top 刪 1/4,508（0.022%）。
10. bottom 刪 21/4,640（0.453%），視覺與 component metadata 均確認為 noise。
11. 沒有可見主體 geometry damage；Local 未被採用。
12. 清理後 56/56 Alpha sweep 有效，且合理鄰域 CV 小幅下降。
13. Alpha Wrap 本身約快 7%，但加入 component 後完整 total 比 voxel-only 高。
14. faces 減少 5.63%。
15. 清理後 sweep 沒有 two-sided wrap；仍保留正式 detector。
16. residual component 前後 volume −0.783 mL；raw voxel 前後 −0.725 mL。
17. matching neighborhood CV 由 0.986% 降為 0.867%，小幅改善。
18. 合理參數鄰域 CV 為 0.867%；固定設定 20-run CV 約 0%。
19. 是；deployment mean/median/P95 = 223.02/210.31/274.27 ms。
20. 最大 bottleneck 是 Alpha Wrap，其次是 radius component search。

### 最終推薦

```text
Recommended Final Pipeline

Voxel:                         1.5 mm
SOR:                           Disabled
ROR:                           Disabled

Connected Component:           Keep dominant main component after dominance/spatial check
Connectivity Radius:           4.0 mm = 2.67 × voxel

Adaptive Local:                Disabled for production
KNN K:                         N/A (diagnostic: density 10, PCA 20)
MAD Factor:                    N/A (diagnostic: 6)

End Protection:                PCA top/bottom 10%; diagnostic local threshold ×3
                                confirmed disconnected clusters are not retained

Original Points:               1,573,158
After Voxel:                   60,129
After SOR:                     60,129 (disabled)
After ROR:                     60,129 (disabled)
After Component:               60,018
Final Points:                  60,018

Total Removed:                 111 component points
Removed Percentage:            0.1846% of voxel cloud

Length Before:                 272.629 mm
Length After:                  270.153 mm
Length Change:                 2.476 mm confirmed noise; 0 mm unexplained damage

Top Removed:                   1 / 4,508 = 0.0222%
Bottom Removed:                21 / 4,640 = 0.4526%

Cleaning Runtime:              component/denoise mean 56.90/59.69 ms
Voxel + cleaning mean:         81.54 ms

Alpha:                         6.25 mm
Offset:                        0.10 mm

Alpha Wrap Runtime:            mean 141.40 ms, P95 166.62 ms
Volume Runtime:                mean 0.014 ms

Deployment Total Runtime:      mean 223.02 ms, median 210.31 ms, P95 274.27 ms

Mesh Vertices:                 1,359
Mesh Faces:                    2,714

Watertight:                    Yes
Manifold:                      Yes
Two-Sided Wrap:                No

Volume:                        388.686 mL
Volume CV:                     0.867% over reasonable parameter neighborhood

Visible Noise Remaining:       No
Geometry Damage:               No evidence
< 500 ms:                      Yes (20/20)
Recommended for real-time cucumber volume: Yes, for memory/stream input
```

分類推薦：Safest、Cleanest、Fastest acceptable、Best Overall 都是 **Voxel + Component**；若允許可見 cluster noise，Voxel-only 才是純計算最快。Adaptive Local 保留為未來資料的診斷功能，不放入本資料 production path。

本輪主要輸出：`results/full_cleaning_alpha_wrap_benchmark.csv`（197 rows、105 columns）、`full_cleaning_repeat_best.csv`、`advanced_*_summary.csv`、`residual_noise_stage_summary.csv`、`advanced_end_cross_section_metrics.csv`、`cleaning/`、`residual_cleaning/`、`meshes/best_*_mesh.{ply,stl}`、`figures/advanced_cleaning/`。
