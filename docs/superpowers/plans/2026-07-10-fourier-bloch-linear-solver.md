# Fourier–Bloch 線形固有値ソルバ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** SNR 線形解析の解像度アーティファクトを、1周期 Fourier 格子＋Bloch シフト（∂_y→D+iK）＋スペクトル補間背景で根本解決する新モジュールを作る。

**Architecture:** 既存 `snr_linear.py`（タイル法・検証済み）は無改変で温存し、新規 `snr_linear_bloch.py` を追加する。摂動を δĝ(y)=e^{iKy}û(y) と分解すると、周期背景では全項の e^{iKy} が消え、演算子は 1周期上で「∂_y→D+iK」の置換だけで得られる。K は明示パラメータとして走査し、背景は平衡の Chebyshev 解から barycentric スペクトル補間で供給する。

**Tech Stack:** Python 3, NumPy, SciPy, pytest, uv（パッケージ管理）。正規化 c=1, ω_pe=1, m_e=1。

## Global Constraints

- パッケージ実行は `uv run ...` を使う。
- 正規化は既存と同一: c=1, ω_pe=1, m_e=1, 長さ c/ω_pe。
- フィールド並びは既存と同一: 各成分 (n,vx,vy) を成分順に並べ、末尾に Ex,Ey,Bz。成分順は `bg.species`（= n の dict 順、SNR では e,inc,ref）。
- 既存 `src/fmi/snr_linear.py` は**改変しない**（回帰基準として温存）。共有物 `_fourier_diff_matrix`,`_field_index`,`BackgroundProfiles` はそこから import して再利用（DRY）。
- 図 PNG は `work/` 直下に置き、docs からは `../work/` 相対参照（コピーしない）。
- コミットは Conventional Commits。デフォルトブランチ上なら作業ブランチを切ってから。

---

### Task 1: Bloch シフト演算子 `build_operator_bloch`

**Files:**
- Create: `src/fmi/snr_linear_bloch.py`
- Test: `tests/test_snr_linear_bloch.py`

**Interfaces:**
- Consumes: `fmi.snr_linear.BackgroundProfiles`, `fmi.snr_linear._fourier_diff_matrix`, `fmi.snr_linear._field_index`
- Produces: `build_operator_bloch(bg: BackgroundProfiles, kx: float, K: float, gamma_ad: float = 5/3, coupling: float = 1.0) -> NDArray` — 複素 (3S+3)M×(3S+3)M 演算子。M=bg.M（1周期）。

- [ ] **Step 1: 失敗するテストを書く（光波 K=0 と Bloch シフト K≠0）**

`tests/test_snr_linear_bloch.py`:

```python
"""Fourier-Bloch 線形固有値ソルバの回帰テスト。"""

import numpy as np
import pytest

from fmi.snr_linear import BackgroundProfiles
from fmi.snr_linear_bloch import build_operator_bloch


def _uniform_single_electron(M: int = 24, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    return BackgroundProfiles(
        y=y, n={"e": np.ones(M)}, beta={"e": 0.0}, T={"e": 0.0},
        mass={"e": 1.0}, charge={"e": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def test_bloch_lightwave_K0():
    """K=0 で Bloch 演算子は光波/Langmuir 分散 ω²=1+kx² を再現。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5, 1.3):
        w2 = np.linalg.eigvals(build_operator_bloch(bg, kx, K=0.0)) ** 2
        assert np.min(np.abs(w2 - (1.0 + kx**2))) < 1e-8
        assert np.min(np.abs(w2 - 1.0)) < 1e-8


def test_bloch_shift_adds_ky():
    """K≠0 は面内波数 ky=K として入り ω²=1+kx²+K² を与える（iK シフトの検証）。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5):
        for K in (0.3, 0.7):
            w2 = np.linalg.eigvals(build_operator_bloch(bg, kx, K=K)) ** 2
            assert np.min(np.abs(w2 - (1.0 + kx**2 + K**2))) < 1e-8
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: FAIL（`ModuleNotFoundError: No module named 'fmi.snr_linear_bloch'`）

- [ ] **Step 3: 最小実装を書く**

`src/fmi/snr_linear_bloch.py`:

```python
"""SNR平衡の線形固有値解析（Fourier-Bloch 版）。

1周期 [0,λ0) 上で摂動を δ(y)=e^{iKy}û(y)（û は λ0 周期）と分解する。周期背景では
全項の e^{iKy} が消え、演算子は「∂_y → D + iK」の置換だけで得られる。K は明示
パラメータ、kx は軸方向波数。タイル化・Bloch事後抽出・線形補間を用いない。

正規化は snr_linear と同一（c=1, ω_pe=1, m_e=1）。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from fmi.snr_linear import BackgroundProfiles, _field_index, _fourier_diff_matrix


def build_operator_bloch(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    gamma_ad: float = 5.0 / 3.0,
    coupling: float = 1.0,
) -> NDArray:
    """ω x = L(kx,K) x の演算子を1周期上で組む（∂_y → D + iK）。"""
    species = bg.species
    M = bg.M
    S = len(species)
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])       # 1周期 = lambda0
    D = _fourier_diff_matrix(M, L_box)
    DK = D + 1j * K * np.eye(M)                   # Bloch シフト
    eye = np.eye(M)
    idx = _field_index(species)
    dim = (3 * S + 3) * M
    Lop = np.zeros((dim, dim), dtype=np.complex128)

    def block(name: str) -> slice:
        b = idx[name]
        return slice(b * M, (b + 1) * M)

    def add(row: str, col: str, mat: NDArray) -> None:
        Lop[block(row), block(col)] += mat

    Ex, Ey, Bz = "Ex", "Ey", "Bz"

    for s in species:
        q = bg.charge[s]
        m = bg.mass[s]
        v0 = bg.beta[s]
        T = bg.T[s]
        n0 = bg.n[s]
        Dn0 = np.diag(n0)
        inv_n0 = np.diag(1.0 / n0)
        DB = np.diag(bg.B0z)
        sn, svx, svy = f"{s}_n", f"{s}_vx", f"{s}_vy"

        # 連続
        add(sn, sn, kx * v0 * eye)
        add(sn, svx, kx * Dn0)
        add(sn, svy, -1j * (DK @ Dn0))

        # 運動 x
        add(svx, svx, kx * v0 * eye)
        add(svx, Ex, 1j * q / m * eye)
        add(svx, svy, 1j * q / m * DB)
        add(svx, sn, (kx * gamma_ad * T / m) * inv_n0)

        # 運動 y
        add(svy, svy, kx * v0 * eye)
        add(svy, Ey, 1j * q / m * eye)
        add(svy, svx, -1j * q / m * DB)
        add(svy, Bz, -1j * q * v0 / m * eye)
        add(svy, sn, -1j * (gamma_ad * T / m) * (inv_n0 @ DK))

        # Ampère（成分電流）
        add(Ex, svx, -1j * coupling * q * Dn0)
        add(Ex, sn, -1j * coupling * q * v0 * eye)
        add(Ey, svy, -1j * coupling * q * Dn0)

    # Faraday
    add(Bz, Ey, kx * eye)
    add(Bz, Ex, 1j * DK)

    # Ampère（∇×B）
    add(Ex, Bz, 1j * DK)
    add(Ey, Bz, kx * eye)

    return Lop


__all__ = ["build_operator_bloch"]
```

- [ ] **Step 4: テストを実行して成功を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: PASS（2 passed）

- [ ] **Step 5: コミット**

```bash
git add src/fmi/snr_linear_bloch.py tests/test_snr_linear_bloch.py
git commit -m "feat(snr_linear_bloch): Bloch-shift operator (d_y -> D+iK) on one period"
```

---

### Task 2: モード抽出 `solve_modes_bloch` と `scan_K_spectrum`

**Files:**
- Modify: `src/fmi/snr_linear_bloch.py`
- Test: `tests/test_snr_linear_bloch.py`

**Interfaces:**
- Consumes: `build_operator_bloch`
- Produces:
  - `solve_modes_bloch(bg, kx, K, n_modes=6, growth_tol=1e-6, **op_kw) -> list[dict]` — 各要素 `{"omega": complex, "gamma": float, "eigvec": NDArray}`（成長率降順、eigvec は û）。
  - `scan_K_spectrum(bg, kx, K_over_k0, growth_tol=1e-7, **op_kw) -> dict` — `{"K_over_k0": NDArray, "gamma": NDArray}`。

- [ ] **Step 1: 失敗するテストを書く（cold filamentation を K=ky として再現）**

`tests/test_snr_linear_bloch.py` に追記:

```python
from fmi.snr_linear_bloch import scan_K_spectrum, solve_modes_bloch


def _two_cold_beams(beta: float, M: int = 48, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    half = 0.5 * np.ones(M)
    return BackgroundProfiles(
        y=y, n={"b1": half, "b2": half}, beta={"b1": beta, "b2": -beta},
        T={"b1": 0.0, "b2": 0.0}, mass={"b1": 1.0, "b2": 1.0},
        charge={"b1": -1.0, "b2": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def _g_formula(ky: float, beta: float) -> float:
    a = ky**2 + 1.0
    return np.sqrt((np.sqrt(a**2 + 4 * beta**2 * ky**2) - a) / 2.0)


def test_bloch_cold_filamentation_via_K():
    """K が面内波数 ky として働き、cold filamentation 成長率 γ(ky=K) を再現。

    一様背景は並進不変なので Bloch 演算子は各 Fourier 高調波 ky=K+n·k0
    (k0=1) に分解され、成長スペクトルは {g_formula(|K+n|)} の集合になる。
    g_formula(K)（最小|ky|）はその集合に「存在する」が最大ではない
    （最大は高|n| の漸近値 β）。よって既存タイル版
    test_cold_filamentation_dispersion と同様に「最大」ではなく「存在」を検査する。
    """
    beta = 0.2
    bg = _two_cold_beams(beta)
    for ky in (0.5, 1.0, 2.0):
        # 全成長モードを集める（g_formula(K) はスペクトル最小側にあるため
        # n_modes は行列サイズ 432 を上回る値にして取りこぼしを防ぐ）
        modes = solve_modes_bloch(bg, kx=0.0, K=ky, n_modes=500, growth_tol=1e-9)
        gammas = np.array([m["gamma"] for m in modes])
        assert gammas.size > 0, f"no growing mode at K={ky}"
        assert np.min(np.abs(gammas - _g_formula(ky, beta))) < 1e-6


def test_scan_K_spectrum_shape():
    """scan_K_spectrum は K_over_k0 と gamma を同長で返す。"""
    beta = 0.2
    bg = _two_cold_beams(beta)
    spec = scan_K_spectrum(bg, kx=0.0, K_over_k0=[0.1, 0.2, 0.3])
    assert spec["K_over_k0"].shape == spec["gamma"].shape == (3,)
    assert np.all(spec["gamma"] >= 0.0)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -k "cold_filamentation_via_K or scan_K_spectrum_shape" -q`
Expected: FAIL（`ImportError: cannot import name 'scan_K_spectrum'`）

- [ ] **Step 3: 最小実装を書く**

`src/fmi/snr_linear_bloch.py` に追記（`__all__` も更新）:

```python
def solve_modes_bloch(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    n_modes: int = 6,
    growth_tol: float = 1e-6,
    **op_kw: float,
) -> list[dict]:
    """L(kx,K) を解き、成長モード（Im ω>growth_tol）を成長率降順で返す。"""
    Lop = build_operator_bloch(bg, kx, K, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    order = np.argsort(-eigval.imag)
    modes: list[dict] = []
    for j in order:
        if eigval[j].imag <= growth_tol:
            break
        modes.append(
            {
                "omega": complex(eigval[j]),
                "gamma": float(eigval[j].imag),
                "eigvec": eigvec[:, j].copy(),
            }
        )
        if len(modes) >= n_modes:
            break
    return modes


def scan_K_spectrum(
    bg: BackgroundProfiles,
    kx: float,
    K_over_k0,
    growth_tol: float = 1e-7,
    **op_kw: float,
) -> dict:
    """Bloch 波数 K を明示走査し、各 K の最大成長率 γ(K) を返す。"""
    k0 = 2.0 * np.pi / bg.lambda0
    r = np.asarray(K_over_k0, dtype=float)
    gamma = np.zeros(r.size)
    for i, ri in enumerate(r):
        modes = solve_modes_bloch(bg, kx, ri * k0, n_modes=1,
                                  growth_tol=growth_tol, **op_kw)
        gamma[i] = modes[0]["gamma"] if modes else 0.0
    return {"K_over_k0": r, "gamma": gamma}
```

`__all__` を更新:

```python
__all__ = ["build_operator_bloch", "solve_modes_bloch", "scan_K_spectrum"]
```

- [ ] **Step 4: テストを実行して成功を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: PASS（4 passed）

- [ ] **Step 5: コミット**

```bash
git add src/fmi/snr_linear_bloch.py tests/test_snr_linear_bloch.py
git commit -m "feat(snr_linear_bloch): solve_modes_bloch and explicit-K scan_K_spectrum"
```

---

### Task 3: スペクトル補間背景 `build_background_1period`

**Files:**
- Modify: `src/fmi/snr_linear_bloch.py`
- Test: `tests/test_snr_linear_bloch.py`

**Interfaces:**
- Consumes: `solve_snr_equilibrium_electron_frame` の返り値 dict、`BackgroundProfiles`
- Produces:
  - `_chebyshev_nodes_weights(N: int) -> tuple[NDArray, NDArray]` — Chebyshev-Lobatto ノード x∈[-1,1]（降順）と barycentric 重み。
  - `_bary_interp(x_nodes, w, f_nodes, x_query) -> NDArray` — Berrut–Trefethen barycentric 補間。
  - `build_background_1period(eq: dict, mime: float, eta: float, points_per_period: int) -> BackgroundProfiles` — 1周期一様格子（n_periods=1）に平衡をスペクトル補間で載せた背景。

- [ ] **Step 1: 失敗するテストを書く（ノード厳密性・スペクトル精度・パリティ）**

`tests/test_snr_linear_bloch.py` に追記:

```python
from fmi.mach_parameters import MachConfig, build_solver_params
from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear_bloch import (
    _bary_interp,
    _chebyshev_nodes_weights,
    build_background_1period,
)


def test_bary_interp_exact_on_polynomial():
    """barycentric 補間は次数<=N の多項式を機械精度で再現。"""
    N = 12
    x, w = _chebyshev_nodes_weights(N)
    f = 3.0 * x**3 - 2.0 * x + 1.0          # 3次（<=N）
    xq = np.linspace(-1.0, 1.0, 37)
    got = _bary_interp(x, w, f, xq)
    exact = 3.0 * xq**3 - 2.0 * xq + 1.0
    assert np.max(np.abs(got - exact)) < 1e-10


def _warm_eq() -> dict:
    return solve_snr_equilibrium_electron_frame(
        eta=0.2, beta_inc=-0.05, beta_ref=0.19,
        Te=0.5, Tinc=0.1, Tref=0.1,
        a0_target=0.25, N_points=64, n_steps=25,
    )


def test_background_1period_basic():
    """1周期背景は n_periods=1・正しい長さ・電荷/質量割当を持つ。"""
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=48)
    assert bg.n_periods == 1
    assert bg.M == 48
    assert set(bg.species) == {"e", "inc", "ref"}
    assert bg.charge == {"e": -1.0, "inc": 1.0, "ref": 1.0}
    assert bg.lambda0 == pytest.approx(float(eq["lambda0"]))
    # B0z は奇対称: y=0 で ~0
    assert abs(bg.B0z[0]) < 1e-3


def test_background_1period_spectral_accuracy():
    """スペクトル補間は解像度を上げても平均密度=1 を保ち滑らか（高調波が小さい）。"""
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=64)
    # 電子密度の空間平均は 1（正規化）
    assert np.mean(bg.n["e"]) == pytest.approx(1.0, abs=1e-3)
    # B0z の最高調波成分は主要成分よりずっと小さい（滑らか=帯域制限）
    sp = np.abs(np.fft.rfft(bg.B0z))
    assert sp[-1] < 1e-3 * sp.max()
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -k "bary_interp or background_1period" -q`
Expected: FAIL（`ImportError: cannot import name '_chebyshev_nodes_weights'`）

- [ ] **Step 3: 最小実装を書く**

`src/fmi/snr_linear_bloch.py` に追記（`__all__` 更新）:

```python
def _chebyshev_nodes_weights(N: int) -> tuple[NDArray, NDArray]:
    """Chebyshev-Gauss-Lobatto ノード x=cos(jπ/N)（降順 1→-1）と barycentric 重み。"""
    x = np.cos(np.arange(N + 1) * np.pi / N)
    w = np.ones(N + 1)
    w[1::2] = -1.0
    w[0] *= 0.5
    w[N] *= 0.5
    return x, w


def _bary_interp(x_nodes: NDArray, w: NDArray, f_nodes: NDArray,
                 x_query) -> NDArray:
    """Berrut–Trefethen barycentric 補間（ノード一致時は値を直接返す）。"""
    xq = np.atleast_1d(np.asarray(x_query, dtype=float))
    out = np.empty(xq.shape, dtype=float)
    for i, xx in enumerate(xq):
        diff = xx - x_nodes
        hit = np.isclose(diff, 0.0)
        if np.any(hit):
            out[i] = f_nodes[int(np.argmax(hit))]
        else:
            t = w / diff
            out[i] = float((t @ f_nodes) / t.sum())
    return out


def build_background_1period(
    eq: dict,
    mime: float,
    eta: float,
    points_per_period: int,
) -> BackgroundProfiles:
    """平衡（Chebyshev 半周期解）を1周期一様格子へ barycentric スペクトル補間で載せる。

    平衡の *_full は半周期 Chebyshev 解の鏡映（B0z は奇、密度は偶）。先頭 N+1 点が
    半周期値。これを Chebyshev barycentric で任意 y に評価し、鏡映パリティを適用して
    [0,λ0) の一様格子を作る（np.interp の C⁰ 折れ点を排除）。
    """
    lambda0 = float(eq["lambda0"])
    L_half = 0.5 * lambda0
    N = (np.asarray(eq["y_full"]).size - 1) // 2
    x_cheb, w_cheb = _chebyshev_nodes_weights(N)

    def half(name: str) -> NDArray:
        return np.asarray(eq[name], dtype=float)[: N + 1]

    ne_h, ni_h, nr_h = half("ne_full"), half("ninc_full"), half("nref_full")
    B_h = half("B0z_full")

    y = np.linspace(0.0, lambda0, points_per_period, endpoint=False)
    s = np.where(y <= L_half, y, lambda0 - y)     # 半周期へ折り返し
    xq = 1.0 - 2.0 * s / L_half                    # y_half=L_half(1-x)/2 の逆写像

    def even(vals_h: NDArray) -> NDArray:
        return _bary_interp(x_cheb, w_cheb, vals_h, xq)

    def odd(vals_h: NDArray) -> NDArray:
        sign = np.where(y <= L_half, 1.0, -1.0)
        return sign * _bary_interp(x_cheb, w_cheb, vals_h, xq)

    params = eq["params"]
    return BackgroundProfiles(
        y=y,
        n={"e": even(ne_h), "inc": even(ni_h), "ref": even(nr_h)},
        beta={
            "e": float(eq["beta_e"]),
            "inc": float(params["beta_inc"]),
            "ref": float(params["beta_ref"]),
        },
        T={
            "e": float(params["Te"]),
            "inc": float(params["Tinc"]),
            "ref": float(params["Tref"]),
        },
        mass={"e": 1.0, "inc": float(mime), "ref": float(mime)},
        charge={"e": -1.0, "inc": 1.0, "ref": 1.0},
        B0z=odd(B_h),
        lambda0=lambda0,
        n_periods=1,
    )
```

`__all__` を更新:

```python
__all__ = [
    "build_operator_bloch",
    "solve_modes_bloch",
    "scan_K_spectrum",
    "build_background_1period",
]
```

- [ ] **Step 4: テストを実行して成功を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: PASS（7 passed）

- [ ] **Step 5: コミット**

```bash
git add src/fmi/snr_linear_bloch.py tests/test_snr_linear_bloch.py
git commit -m "feat(snr_linear_bloch): spectral (barycentric) one-period background"
```

---

### Task 4: タイル法との等価性（warm 平衡）

**Files:**
- Test: `tests/test_snr_linear_bloch.py`

**Interfaces:**
- Consumes: `fmi.snr_linear.build_background`, `fmi.snr_linear.scan_bloch_spectrum`, `fmi.snr_linear.solve_modes`, `build_background_1period`, `scan_K_spectrum`, `solve_modes_bloch`

- [ ] **Step 1: 失敗するテストを書く（warm で γ が一致）**

`tests/test_snr_linear_bloch.py` に追記:

```python
from fmi.snr_linear import build_background as build_background_tiled
from fmi.snr_linear import scan_bloch_spectrum, solve_modes


def test_equivalence_tiled_vs_bloch_warm_fmi():
    """warm 平衡・kx=0 の FMI 最大成長率がタイル法と1周期Bloch法で一致。"""
    eq = _warm_eq()
    # タイル法（K は事後抽出）: FMI バンドの最大成長率
    bg_t = build_background_tiled(eq, mime=400.0, eta=0.2,
                                  points_per_period=32, n_periods=6)
    g_tiled = scan_bloch_spectrum(bg_t, kx=0.0)["gamma"].max()
    # 1周期Bloch法: K を明示走査した最大成長率
    bg_b = build_background_1period(eq, mime=400.0, eta=0.2,
                                    points_per_period=32)
    g_bloch = scan_K_spectrum(
        bg_b, kx=0.0, K_over_k0=np.linspace(0.0, 0.5, 7)
    )["gamma"].max()
    assert abs(g_tiled - g_bloch) / abs(g_tiled) < 1e-2


def test_equivalence_tiled_vs_bloch_warm_oblique():
    """warm 平衡・斜め kx/k0=1.10 の最大成長率が両手法で一致。"""
    eq = _warm_eq()
    k0 = 2.0 * np.pi / float(eq["lambda0"])
    bg_t = build_background_tiled(eq, mime=400.0, eta=0.2,
                                  points_per_period=32, n_periods=4)
    g_tiled = solve_modes(bg_t, kx=1.10 * k0, n_modes=1)[0]["gamma"]
    bg_b = build_background_1period(eq, mime=400.0, eta=0.2,
                                    points_per_period=32)
    g_bloch = solve_modes_bloch(bg_b, kx=1.10 * k0, K=0.0, n_modes=1)[0]["gamma"]
    assert abs(g_tiled - g_bloch) / abs(g_tiled) < 1e-2
```

- [ ] **Step 2: テストを実行して失敗するか確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -k equivalence -q`
Expected: 実装が正しければ PASS。もし FAIL する場合は許容差ではなく実装バグの兆候（背景サンプリングや K 走査範囲）。まず PASS を確認する。

補足: この2テストは新規実装ではなく既存機能の突き合わせなので、Step 2 で PASS すれば Step 3-4 は不要。FAIL したら原因（例: 斜めモードの最大が K≠0 セクタにある）を調べ、`solve_modes_bloch` の K を FMI バンド最速 K に合わせる等でテスト条件を厳密化する。

- [ ] **Step 3: コミット**

```bash
git add tests/test_snr_linear_bloch.py
git commit -m "test(snr_linear_bloch): warm-regime equivalence with tiled Fourier method"
```

---

### Task 5: 固有関数の再構成 `reconstruct_eigenfunction`

**Files:**
- Modify: `src/fmi/snr_linear_bloch.py`
- Test: `tests/test_snr_linear_bloch.py`

**Interfaces:**
- Consumes: `solve_modes_bloch`, `_field_index`
- Produces: `reconstruct_eigenfunction(mode: dict, bg: BackgroundProfiles, K: float, field: str = "Bz", n_display: int = 6) -> tuple[NDArray, NDArray]` — 表示用に δ(y)=e^{iKy}û(y) を n_display 周期展開した (ys, delta)。

- [ ] **Step 1: 失敗するテストを書く（Bloch 位相性）**

`tests/test_snr_linear_bloch.py` に追記:

```python
from fmi.snr_linear_bloch import reconstruct_eigenfunction


def test_reconstruct_bloch_phase():
    """再構成した δ は Bloch 性 δ(y+λ0)=e^{iKλ0}δ(y) を満たす。"""
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=32)
    k0 = 2.0 * np.pi / bg.lambda0
    K = 0.3 * k0
    mode = solve_modes_bloch(bg, kx=0.0, K=K, n_modes=1)[0]
    ys, delta = reconstruct_eigenfunction(mode, bg, K, field="Bz", n_display=3)
    M = bg.M
    assert ys.size == delta.size == 3 * M
    # 1周期ずらすと位相因子 e^{iKλ0} 倍
    phase = np.exp(1j * K * bg.lambda0)
    assert np.allclose(delta[M:2 * M], phase * delta[0:M], atol=1e-10)
```

- [ ] **Step 2: テストを実行して失敗を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -k reconstruct -q`
Expected: FAIL（`ImportError: cannot import name 'reconstruct_eigenfunction'`）

- [ ] **Step 3: 最小実装を書く**

`src/fmi/snr_linear_bloch.py` に追記（`__all__` 更新）:

```python
def reconstruct_eigenfunction(
    mode: dict,
    bg: BackgroundProfiles,
    K: float,
    field: str = "Bz",
    n_display: int = 6,
) -> tuple[NDArray, NDArray]:
    """固有ベクトル û から δ(y)=e^{iKy}û(y) を n_display 周期展開して返す。"""
    M = bg.M
    lam = bg.lambda0
    b = _field_index(bg.species)[field]
    uhat = mode["eigvec"][b * M:(b + 1) * M]
    ys = np.concatenate([bg.y + k * lam for k in range(n_display)])
    delta = np.exp(1j * K * ys) * np.tile(uhat, n_display)
    return ys, delta
```

`__all__` に `"reconstruct_eigenfunction"` を追加。

- [ ] **Step 4: テストを実行して成功を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: PASS（10 passed）

- [ ] **Step 5: コミット**

```bash
git add src/fmi/snr_linear_bloch.py tests/test_snr_linear_bloch.py
git commit -m "feat(snr_linear_bloch): reconstruct_eigenfunction (delta = e^{iKy} uhat)"
```

---

### Task 6: warm 域での平滑固有関数の検証とデモ図

**設計変更（2026-07-11）:** 当初は cold σ=0 で「アーティファクト消失」を示す予定だったが、実測で cold の最速成長モードは物理 FMI ではなく **spurious（純 ref_vy、δn≈δE≈δB≈機械精度ゼロ、n_ref≈5.9e-54 のセルに局在）** と判明。原因は運動yの `inv_n0=1/n0` が cold で密度→0 の領域で発散すること（タイル法も同じ挙動＝Bloch のバグではない）。この inv_n0 縮退は本分岐の対象外（別途 spec が必要）。よって Task 6 は **warm 平衡**（密度最小 0.12、inv_n0 健全）で「Bloch 法が物理 FMI モードの滑らかな δB_z 固有関数を与える」ことを検証する。cold の限界は docs に既知限界として明記する。

**Files:**
- Test: `tests/test_snr_linear_bloch.py`
- Create: `work/plot_snr_bloch_demo.py`

**Interfaces:**
- Consumes: すべての公開 API、`_warm_eq`（既存ヘルパ）、`_field_index`

**実測（warm, kx=0, K=0 の FMI モード）:** γ=5.843e-4（PPP=24/32/48/64 で rel 7.7e-14 収束）。場構成は **Bz 支配（0.662）**, ref_n(0.515), inc_n(0.51), ref_vy(0.168) の磁気 filamentation モード。δB_z の Nyquist/最大成分 = 2.4e-13（帯域制限＝滑らか）。

- [ ] **Step 1: 失敗するテストを書く（warm FMI の収束と平滑性）**

`tests/test_snr_linear_bloch.py` に追記（`_field_index` をファイル冒頭の `fmi.snr_linear` import に追加）:

```python
def test_warm_bloch_gamma_converges():
    """warm FMI（kx=0,K=0）成長率が PPP=32→64 で収束（平滑背景＋スペクトル補間）。"""
    eq = _warm_eq()

    def gamma_at(ppp: int) -> float:
        bg = build_background_1period(eq, mime=400.0, eta=0.2,
                                      points_per_period=ppp)
        m = solve_modes_bloch(bg, kx=0.0, K=0.0, n_modes=1, growth_tol=1e-7)
        return m[0]["gamma"] if m else 0.0

    g32, g64 = gamma_at(32), gamma_at(64)
    assert g32 > 1e-4 and g64 > 1e-4
    assert abs(g32 - g64) / abs(g64) < 1e-2      # 実測 ~1e-13


def test_warm_bloch_eigenfunction_smooth():
    """warm FMI（kx=0,K=0）は δB_z 支配の磁気モードで、固有関数が滑らか（帯域制限）。"""
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=48)
    m = solve_modes_bloch(bg, kx=0.0, K=0.0, n_modes=1, growth_tol=1e-7)[0]
    M = bg.M
    v = m["eigvec"]

    def blk_norm(field: str) -> float:
        i = _field_index(bg.species)[field]
        return float(np.linalg.norm(v[i * M:(i + 1) * M]))

    # δB_z が支配的＝磁気 filamentation モード（spurious な純 ref_vy でない）
    assert blk_norm("Bz") == max(blk_norm(f) for f in _field_index(bg.species))
    # 複素固有ベクトルは fft（rfft は実数専用）。Nyquist 成分が主要成分よりずっと小
    bz = v[_field_index(bg.species)["Bz"] * M:(_field_index(bg.species)["Bz"] + 1) * M]
    sp = np.abs(np.fft.fft(bz))
    assert sp[M // 2] < 1e-6 * sp.max()          # Gibbs リンギングなし（実測 2.4e-13）
```

- [ ] **Step 2: テストを実行して確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -k warm_bloch -q`
Expected: PASS（既存 API のみ使用）。

- [ ] **Step 3: デモ図スクリプトを作る**

`work/plot_snr_bloch_demo.py`:

```python
"""Fourier-Bloch ソルバのデモ: warm 平衡の物理 FMI モードの δB_z 固有関数が
滑らか（帯域制限・Gibbs 跳びなし）であることを示し、明示 K 走査の FMI バンドを併記する。

出力は work/ 直下（3s_cf_SNR/FMI の慣習）。cold σ=0 は inv_n0 縮退で最速モードが
spurious になるため、物理モードが健全な warm 平衡を用いる。
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear_bloch import (
    build_background_1period,
    reconstruct_eigenfunction,
    scan_K_spectrum,
    solve_modes_bloch,
)


def main() -> None:
    eq = solve_snr_equilibrium_electron_frame(
        eta=0.2, beta_inc=-0.05, beta_ref=0.19,
        Te=0.5, Tinc=0.1, Tref=0.1,
        a0_target=0.25, N_points=64, n_steps=25,
    )
    lam = float(eq["lambda0"])

    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=48)
    mode = solve_modes_bloch(bg, kx=0.0, K=0.0, n_modes=1, growth_tol=1e-7)[0]
    ys, dBz = reconstruct_eigenfunction(mode, bg, K=0.0, field="Bz", n_display=6)
    peak = np.abs(dBz).max()
    dBz = dBz / peak
    spec = scan_K_spectrum(bg, kx=0.0, K_over_k0=np.linspace(0.0, 0.5, 11))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 3.6))
    axL.plot(ys / lam, dBz.real, color="#c0392b", lw=1.3, label="Re δB_z")
    axL.plot(ys / lam, dBz.imag, color="#c0392b", lw=0.9, ls="--", label="Im δB_z")
    axL.set_xlabel("y/λ0"); axL.set_ylabel("δB_z (norm., Bloch spectral bg)")
    axL.set_title(f"warm FMI eigenfunction (γ={mode['gamma']:.3e})")
    axL.legend(fontsize=8)
    axR.plot(spec["K_over_k0"], spec["gamma"], "o-", color="#2980b9", ms=4)
    axR.set_xlabel("K/k0"); axR.set_ylabel("γ/ω_pe")
    axR.set_title("FMI band (explicit K scan)")
    fig.tight_layout()
    fig.savefig("work/snr_bloch_demo.png", dpi=200)
    print("saved: work/snr_bloch_demo.png")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: デモを実行して図生成を確認**

Run: `uv run python work/plot_snr_bloch_demo.py`
Expected: `saved: work/snr_bloch_demo.png`（エラーなく完了）

- [ ] **Step 5: 全テストを実行して緑を確認**

Run: `uv run pytest tests/test_snr_linear_bloch.py -q`
Expected: PASS（12 passed）

- [ ] **Step 6: コミット**

```bash
git add tests/test_snr_linear_bloch.py work/plot_snr_bloch_demo.py
git commit -m "test(snr_linear_bloch): warm-regime smooth eigenfunction + demo figure"
```

---

## 完了後

- 既存 `tests/test_snr_linear.py`（タイル法の回帰）が引き続き緑であることを確認: `uv run pytest tests/test_snr_linear.py -q`（5/5 pass, 2026-07-11）。
- `docs/snr_linear_analysis.html` に Bloch 法の節を追加（別作業）。warm の平滑 δB_z 固有関数図 `../work/snr_bloch_demo.png` を掲載。
- durable memory `snr-current-filament-linear-analysis` に「Fourier-Bloch 版：warm で物理FMIの平滑固有関数を再現・タイル法と<1%一致」を追記（別作業）。

## 既知の限界（本分岐の対象外）

- **cold σ=0 の inv_n0 縮退**: cold（T_ion→0）ではイオン密度が n_ref≈5.9e-54 まで evacuate し、運動yの圧力項 `inv_n0=1/n0` が発散する。その結果、最速成長モードは物理 FMI ではなく **spurious な純 ref_vy モード**（δn≈δE≈δB≈機械精度ゼロ、固有ペア残差~1e-15 の真の固有値）になる。タイル法 `snr_linear` も同一挙動＝Bloch 法固有のバグではない。Fourier-Bloch は補間 Gibbs は除去するが、この密度→0 の 1/n₀ 縮退は除去しない。将来対応（別 spec）: (a) inv_n0 の正則化（n0 フロア／密度重み弱形式）、(b) spurious フィルタ（|δn|/|δv| 閾値）、(c) cold 非対称平衡の収束性検証（native 残差・Harris 比較）。
