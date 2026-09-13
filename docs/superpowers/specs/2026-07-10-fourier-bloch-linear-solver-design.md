# Fourier–Bloch 線形固有値ソルバ（II′）設計

**日付:** 2026-07-10
**対象:** SNR 電流フィラメント平衡の線形固有値解析における解像度アーティファクトの根本解決

## 背景と動機

現行の線形解析 `src/fmi/snr_linear.py` は、Chebyshev＋鏡映で解いた平衡を **一様 Fourier 格子へ `np.interp`（区分線形 = C⁰）で載せ替え、N 周期にタイル化**して大領域固有値解析を行う。この「2つの異なる離散化の橋渡し」が cold（M_S≈155、薄い電流シート）で固有関数に **Gibbs 振動（見かけの不連続）** を生み、成長率 γ も points-per-period（PPP）24→48 で 2–5% ドリフトする。

実測（2026-07-10）:
- cold: シートは y/λ₀ = 0, 0.5, 1.0 の3面。一様 Fourier はこれを解像しきれない。γ は PPP=16/24/48/72/96 で 6.85→6.60e-3（約4%）ドリフト、収束は PPP≳48。
- warm「natural」(Te=0.5,Ti=0.1): 平衡が滑らかで γ は PPP=24/32/48 で5桁一致（アーティファクトなし）。

**真因:** cold そのものではなく、平衡（Chebyshev）と線形（Fourier）の**離散化不整合**（補間の C⁰ 折れ点＋タイル化）。

**方針転換の指摘（研究指導）:** 最初から周期スペクトル（Fourier）で通せば周期境界条件が自然に入る。→ 本設計は線形段を Fourier-native に統一し、Bloch 位相をシフトで取り込む **II′** を採る。

## 目標 / 非目標

**目標:**
- 線形段を1つの周期 Fourier 枠に統一し、補間・タイル化を撤廃する。
- 固有関数を band-limited（原理的に Gibbs なし）にし、cold の「不連続」を消す。
- Bloch 波数 K を明示パラメータとして走査する（事後抽出を廃止）。
- warm 域で現行 `snr_linear` とスペクトル精度で一致し、cold で γ が M 単調収束することを示す。

**非目標（YAGNI / 将来）:**
- 周期座標マップ（cold の DOF 効率化）。まず素の Fourier-Bloch で跳び消失を確認してから検討。
- 平衡ソルバの Fourier 書き直し（案 II）。
- 既存 `snr_linear.py` は無改変で回帰基準として温存。

## 中核となる定式化：Bloch シフト

1周期 [0,λ₀) 上で摂動を

```
δĝ(y) = e^{iKy} û(y),      û(y+λ₀) = û(y)   （λ₀ 周期）
```

と分解する（δĝ(y) は既存の摂動振幅 δg(x,y,t)=δĝ(y) exp(i kx x − iωt) の y 部分）。周期背景に対し

```
∂_y δĝ = e^{iKy} (∂_y + iK) û
```

全項に共通の e^{iKy} が括り出せて消えるので、**演算子は「∂_y → D + iK·I」の置換だけ**で得られる。ここで D は1周期 [0,λ₀) 上の Fourier 微分行列。

- K はパラメータ（走査）。kx（軸方向）は従来通り別パラメータ。演算子は L(kx, K)。
- 行列サイズは (3S+3)·M（M = points-per-period、1周期分）。従来のタイル ×N より N 倍小さい。
- Bloch 事後抽出（`bloch_wavenumber`）は不要（K は入力）。

**置換対象は D が現れる箇所のみ:** 連続式の ∂_y(n₀δv)、運動 y の圧力 ∂_y(δp)、Faraday の ∂_yδE_x、Ampère の ∇×B。代数項（kx·v0、背景の対角乗算など）は不変。参照: `src/fmi/snr_linear.py::build_operator`。

## 背景の供給（不整合の除去）

`np.interp`（C⁰）を廃止し、平衡の Chebyshev 解を **barycentric スペクトル補間**で1周期一様格子 M 点へ載せる（C^∞）。平衡は滑らか・周期的なので、これで背景の折れ点が消え、Fourier 演算子が滑らかな周期関数を見る。結果、固有関数は band-limited となり解像度不足でも滑らか（跳びは出ない）。

実装メモ: 平衡は半周期 Chebyshev（`x_std=cos(jπ/N)`、`y_half=L_half(1−x_std)/2`）の解 A(x), Φ(x) と、その鏡映（B,E は奇、密度は偶）で全周期を成す。1周期一様 y へは、半周期 Chebyshev 上の barycentric 補間で任意 y を評価し、鏡映パリティを適用して [0,λ₀) を構成する（`np.interp` を barycentric に置換するのが最小差分）。

## コンポーネント / ファイル

新規 `src/fmi/snr_linear_bloch.py`（既存 `snr_linear.py` は温存）:

- `build_background_1period(eq, mime, eta, points_per_period) -> BackgroundProfiles`
  - 既存 `BackgroundProfiles`（`n_periods=1`）を再利用。背景はスペクトル補間で供給。
- `build_operator_bloch(bg, kx, K, gamma_ad=5/3, coupling=1.0) -> NDArray`
  - `build_operator` の写しに D→(D+iK·I) を適用。
- `solve_modes_bloch(bg, kx, K, n_modes=6, growth_tol=1e-6, **op_kw) -> list[dict]`
  - eig → Im ω>tol を降順。返り値 `{"omega","gamma","eigvec"}`（eigvec は û）。
- `scan_K_spectrum(bg, kx, K_values, **op_kw) -> dict`
  - K 明示走査。出力キーは現行互換 `{"K_over_k0","gamma"}`。
- `reconstruct_eigenfunction(mode, bg, K, n_display=6) -> tuple[NDArray, NDArray]`
  - δ(y)=e^{iKy}û(y) を n 周期展開して (y, δ) を返す（表示用の厳密再構成）。

共有ユーティリティ `_fourier_diff_matrix(M, L)` は `snr_linear` から再利用（またはコピー）。フィールド並び `[e_n,e_vx,e_vy, inc_n,inc_vx,inc_vy, ref_n,ref_vx,ref_vy, Ex,Ey,Bz]` は現行と同一。

## データフロー

```
eq(Chebyshev) ─spectral interp─▶ 1周期 uniform M点の背景 (build_background_1period)
             ─build_operator_bloch(kx,K)─▶ ω x = L(kx,K) x（密 eig）
             ─▶ γ(kx,K), û ─reconstruct_eigenfunction─▶ δ(y)=e^{iKy}û（表示）
```

## 検証（TDD、平衡と独立）

新規 `tests/test_snr_linear_bloch.py`:

1. **光波 / Langmuir**: 一様単一電子・K=0・kx 掃引 → ω²=1+kx²（機械精度）。
2. **cold filamentation**: 一様逆流2 cold ビーム・kx=0・K 掃引（K=ky） → 解析 γ²(ky)=½[√((ky²+1)²+4β²ky²)−(ky²+1)] に一致、ky→∞ で γ→β。
3. **等価性（本質テスト）**: warm 平衡で現行 `snr_linear`（タイル）と `snr_linear_bloch`（1周期）の γ(K)・γ_max(kx) がスペクトル精度で一致。
4. **cold アーティファクト消失**: cold の固有関数が滑らか（δB_z の高周波リンギングが閾値以下）、γ が M 増加で単調収束（現行の 2–5% ふらつきが消える）。

## 検証コマンド

```bash
cd /Users/koki/3s_cf_SNR/FMI
uv run pytest tests/test_snr_linear_bloch.py -q
uv run python work/plot_snr_bloch_demo.py   # cold/warm の固有関数と γ(K) を可視化
```

## リスクと対応

- **cold 薄シートの解像度**: 素の一様 Fourier は cold で多モードを要する。ただし出力は「滑らかだが未解像」で Gibbs の跳びではない。真の cold 精度が要る場合は将来の周期座標マップで対応（本設計外）。
- **barycentric 補間の実装**: 鏡映グリッド上の評価はパリティ適用で処理。半周期 Chebyshev 補間器を明示的に持つのが安全。
- **等価性テストの許容差**: タイル法と1周期法は数学的に等価だが有限解像度で微差。warm（滑らか）で比較し相対差 <1e-3 を基準にする。

## 段階的実装方針

1. `_fourier_diff_matrix` と `build_operator_bloch`（D→D+iK）を先に作り、テスト1（光波、K=0）で較正。
2. テスト2（cold filamentation、K=ky）で K シフトの物理を確認。
3. `build_background_1period`（barycentric）を実装、テスト3（warm 等価性）。
4. `reconstruct_eigenfunction` と可視化、テスト4（cold 跳び消失）。
