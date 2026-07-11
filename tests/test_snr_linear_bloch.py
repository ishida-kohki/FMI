"""Fourier-Bloch 線形固有値ソルバの回帰テスト。"""

import numpy as np
import pytest

from fmi.mach_parameters import MachConfig, build_solver_params
from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear import BackgroundProfiles
from fmi.snr_linear_bloch import (
    _bary_interp,
    _chebyshev_nodes_weights,
    build_background_1period,
    build_operator_bloch,
    scan_K_spectrum,
    solve_modes_bloch,
)


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
    g_formula(K)（最小|ky|）はその集合に存在するが最大ではない（最大は
    高|n| の漸近値 β）。既存タイル版 test_cold_filamentation_dispersion と
    同様に「最大」ではなく「存在」を検査する。
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
