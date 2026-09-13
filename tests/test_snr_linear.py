"""SNR線形固有値解析の回帰テスト。

演算子を平衡と独立な2つの解析解（光波分散・cold filamentation）で較正し、
SNR平衡背景で FMI 不安定帯（Bloch ゾーン全域で成長）と解像度収束を確認する。
"""

import numpy as np
import pytest

from fmi.mach_parameters import MachConfig, build_solver_params
from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear import (
    BackgroundProfiles,
    _field_index,
    _fourier_diff_matrix,
    build_background,
    build_operator,
    scan_bloch_spectrum,
    solve_modes,
)


def _uniform_single_electron(M: int = 24, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    return BackgroundProfiles(
        y=y, n={"e": np.ones(M)}, beta={"e": 0.0}, T={"e": 0.0},
        mass={"e": 1.0}, charge={"e": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def _two_cold_beams(beta: float, M: int = 48, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    half = 0.5 * np.ones(M)
    return BackgroundProfiles(
        y=y, n={"b1": half, "b2": half}, beta={"b1": beta, "b2": -beta},
        T={"b1": 0.0, "b2": 0.0}, mass={"b1": 1.0, "b2": 1.0},
        charge={"b1": -1.0, "b2": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def _snr_background(points_per_period: int, n_periods: int) -> BackgroundProfiles:
    # 本プロジェクト既定の非磁化（σ=0、音速マッハ基準）平衡を背景に使う。
    sp = build_solver_params(
        MachConfig(M_S=154.92, beta_sh=0.25, beta_te=0.1, sigma=0.0)
    )
    eq = solve_snr_equilibrium_electron_frame(
        **sp, a0_target=0.25, N_points=64, n_steps=25
    )
    return build_background(
        eq, mime=400.0, eta=0.2,
        points_per_period=points_per_period, n_periods=n_periods,
    )


def test_lightwave_calibration():
    """一様単一電子で光波/Langmuir 分散 ω²=ω_pe²+c²kx² (=1+kx²) を再現。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5, 1.3):
        w2 = np.linalg.eigvals(build_operator(bg, kx)) ** 2
        assert np.min(np.abs(w2 - (1.0 + kx**2))) < 1e-8   # k_y=0 横波
        assert np.min(np.abs(w2 - 1.0)) < 1e-8             # Langmuir


def test_cold_filamentation_dispersion():
    """逆流2ビームの cold filamentation 成長率が解析分散と一致。"""
    beta = 0.2
    bg = _two_cold_beams(beta)          # L=2π → k_y=整数
    w = np.linalg.eigvals(build_operator(bg, kx=0.0))
    gammas = w.imag[w.imag > 1e-9]

    def g_formula(ky: float) -> float:
        a = ky**2 + 1.0
        return np.sqrt((np.sqrt(a**2 + 4 * beta**2 * ky**2) - a) / 2.0)

    for ky in (1.0, 2.0, 3.0):
        gf = g_formula(ky)
        assert np.min(np.abs(gammas - gf)) < 1e-6   # 予測成長率が固有値に存在
    # k_y→∞ の漸近: 最大成長率 → β ω_pe
    assert gammas.max() == pytest.approx(beta, rel=0.05)


# 低温（密度空乏あり）背景の線形解析には density_floor>0 が実質必須:
# 運動 y の 1/n0（圧力項・背景場力項）が空乏域で発散するため。
_COLD_FLOOR = 3e-2


def test_snr_fmi_band_unstable():
    """SNR平衡・kx=0 で Bloch ゾーン全域が不安定、併合枝 K=k0/2 も成長。"""
    bg = _snr_background(points_per_period=20, n_periods=8)
    spec = scan_bloch_spectrum(bg, kx=0.0, density_floor=_COLD_FLOOR)
    K = spec["K_over_k0"]
    g = spec["gamma"]
    assert g.size > 0
    assert np.all(g > 1e-3)                       # 全セクタ不安定（広帯域）
    # 併合（period-doubling）セクタ K/k0=0.5 が解像され不安定
    i_merge = int(np.argmin(np.abs(K - 0.5)))
    assert abs(K[i_merge] - 0.5) < 1e-6
    assert g[i_merge] > 1e-3


def test_snr_oblique_growth_rate_converges():
    """斜め枝の最大成長率が解像度(points_per_period)に対し収束（spurious でない）。

    cold は PPP≳36 が必要（実測: PPP36 9.21e-3, PPP48 9.12e-3, 差 0.9%）。
    """
    def g_at(ppp: int) -> float:
        bg = _snr_background(points_per_period=ppp, n_periods=4)
        k0 = 2.0 * np.pi / bg.lambda0
        return solve_modes(bg, kx=0.3 * k0, n_modes=1,
                           density_floor=_COLD_FLOOR)[0]["gamma"]

    g1, g2 = g_at(36), g_at(48)
    assert abs(g1 - g2) / abs(g2) < 0.05


def test_snr_oblique_mode_faster_than_fmi():
    """最速モードは有限 kx（斜め/ドリフトキンク枝）で kx=0 の FMI より明確に速い。

    補正後（背景場力項込み）の実測: γ_obl≈9.2e-3 vs γ_FMI≈2.5e-3（PPP=36, floor=3e-2）。
    """
    bg = _snr_background(points_per_period=36, n_periods=4)
    k0 = 2.0 * np.pi / bg.lambda0
    g0 = solve_modes(bg, kx=0.0, n_modes=1,
                     density_floor=_COLD_FLOOR)[0]["gamma"]
    g_obl = solve_modes(bg, kx=0.3 * k0, n_modes=1,
                        density_floor=_COLD_FLOOR)[0]["gamma"]
    assert g0 > 1e-3                       # kx=0 の FMI も不安定
    assert g_obl > 1.3 * g0                # 斜め枝が明確に速い（実測 ~3.6x）


def test_cold_floor_stable_magnetic_mode():
    """低温斜め枝は floor に頑健で、磁場に結合した物理モード。

    密度下限（floor）を 2 桁振っても γ の変動が小さく、固有ベクトルは
    δB_z に有意に結合する（場と非結合な spurious 表現でない）。
    """
    bg = _snr_background(points_per_period=24, n_periods=4)
    k0 = 2.0 * np.pi / bg.lambda0
    idx = _field_index(bg.species)
    M = bg.M

    gammas = []
    bz_fracs = []
    for fl in (1e-3, 1e-2, 3e-2, 1e-1):
        m = solve_modes(bg, kx=0.3 * k0, n_modes=1, density_floor=fl)[0]
        gammas.append(m["gamma"])
        norms = {f: float(np.linalg.norm(m["eigvec"][i * M:(i + 1) * M]))
                 for f, i in idx.items()}
        bz_fracs.append(norms["Bz"] / max(norms.values()))
    g = np.asarray(gammas)
    assert (g.max() - g.min()) / g.min() < 0.05   # floor 2桁で γ 変動 <5%
    assert min(bz_fracs) > 0.1                    # δB_z に有意に結合


def test_density_floor_inactive_when_resolved():
    """密度が空乏化しない warm 平衡では floor が非作動（γ 不変）。"""
    sp = dict(eta=0.2, beta_inc=-0.05, beta_ref=0.19,
              Te=0.5, Tinc=0.1, Tref=0.1)
    eq = solve_snr_equilibrium_electron_frame(
        **sp, a0_target=0.25, N_points=64, n_steps=25
    )
    bg = build_background(eq, mime=400.0, eta=0.2,
                          points_per_period=24, n_periods=4)
    k0 = 2.0 * np.pi / bg.lambda0
    g_raw = solve_modes(bg, kx=1.10 * k0, n_modes=1)[0]["gamma"]
    g_reg = solve_modes(bg, kx=1.10 * k0, n_modes=1,
                        density_floor=3e-2)[0]["gamma"]
    assert abs(g_reg - g_raw) / g_raw < 1e-6


def test_translation_zero_mode_annihilated():
    """平衡の y 微小並進（ω=0 の厳密解）を演算子が消す（背景場力項の検査）。

    並進モード δn_s=-n0_s', δEy=-E0y', δBz=-B0z', δv=0 は kx=0 で L x = 0 を
    満たすべき。運動 y の背景場力項 q_s δn_s (E0y - β_s B0z) が欠けていると
    運動 y ブロックに O(1e-4) の残差が残る（欠落時実測 9.2e-5, 補正後 4.6e-6）。
    γ_ad=1（等温、閉包が背景 Boltzmann と自己無撞着）で検査する。
    """
    sp = dict(eta=0.2, beta_inc=-0.05, beta_ref=0.19,
              Te=0.5, Tinc=0.1, Tref=0.1)
    eq = solve_snr_equilibrium_electron_frame(
        **sp, a0_target=0.25, N_points=64, n_steps=25
    )
    bg = build_background(eq, mime=400.0, eta=0.2,
                          points_per_period=32, n_periods=2)
    M = bg.M
    S = len(bg.species)
    idx = _field_index(bg.species)
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])
    D = _fourier_diff_matrix(M, L_box)

    x = np.zeros((3 * S + 3) * M, dtype=np.complex128)
    for s in bg.species:
        x[idx[f"{s}_n"] * M:(idx[f"{s}_n"] + 1) * M] = -(D @ bg.n[s])
    x[idx["Ey"] * M:(idx["Ey"] + 1) * M] = -(D @ bg.E0y)
    x[idx["Bz"] * M:(idx["Bz"] + 1) * M] = -(D @ bg.B0z)

    r = build_operator(bg, kx=0.0, gamma_ad=1.0) @ x
    xn = np.linalg.norm(x)
    for f in ("e_vy", "inc_vy", "ref_vy"):
        i = idx[f]
        assert np.linalg.norm(r[i * M:(i + 1) * M]) / xn < 1e-5
