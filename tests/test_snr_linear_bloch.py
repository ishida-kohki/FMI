"""Fourier-Bloch 線形固有値ソルバの回帰テスト。"""

import numpy as np
import pytest

from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear import BackgroundProfiles, _field_index, scan_bloch_spectrum, solve_modes
from fmi.snr_linear import build_background as build_background_tiled
from fmi.snr_linear_bloch import (
    _bary_interp,
    _chebyshev_nodes_weights,
    build_background_1period,
    build_background_fourier,
    build_operator_bloch,
    reconstruct_eigenfunction,
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


def test_warm_fmi_marginal_both_methods():
    """warm・kx=0 の FMI は補正後ほぼ安定化し、両手法とも γ<5e-5 で一致。

    背景場力項 q_s δn_s (E0y - β_s B0z) の追加で warm FMI（補正前 γ≈5.8e-4）は
    実質安定化する（tiled 実測 ~1e-5、bloch は成長モードなし）。この周辺性を
    両手法でロックする。
    """
    eq = _warm_eq()
    bg_t = build_background_tiled(eq, mime=400.0, eta=0.2,
                                  points_per_period=32, n_periods=6)
    g_t = scan_bloch_spectrum(bg_t, kx=0.0)["gamma"]
    g_tiled = float(g_t.max()) if g_t.size else 0.0
    bg_b = build_background_1period(eq, mime=400.0, eta=0.2,
                                    points_per_period=32)
    g_bloch = scan_K_spectrum(
        bg_b, kx=0.0, K_over_k0=np.linspace(0.0, 0.5, 7)
    )["gamma"].max()
    assert g_tiled < 5e-5
    assert g_bloch < 5e-5


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


def test_reconstruct_bloch_phase():
    """再構成した δ は Bloch 性 δ(y+λ0)=e^{iKλ0}δ(y) を満たす。

    kx=0 の FMI は補正後ほぼ安定なので、成長モードが確実に存在する
    斜め枝（kx/k0=1.10）で検査する。
    """
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=32)
    k0 = 2.0 * np.pi / bg.lambda0
    K = 0.3 * k0
    mode = solve_modes_bloch(bg, kx=1.10 * k0, K=K, n_modes=1)[0]
    ys, delta = reconstruct_eigenfunction(mode, bg, K, field="Bz", n_display=3)
    M = bg.M
    assert ys.size == delta.size == 3 * M
    # 1周期ずらすと位相因子 e^{iKλ0} 倍
    phase = np.exp(1j * K * bg.lambda0)
    assert np.allclose(delta[M:2 * M], phase * delta[0:M], atol=1e-10)


def test_warm_bloch_gamma_converges():
    """warm 斜め枝（kx/k0=1.10, K=0）成長率が PPP=32→64 で収束。"""
    eq = _warm_eq()
    k0 = 2.0 * np.pi / float(eq["lambda0"])

    def gamma_at(ppp: int) -> float:
        bg = build_background_1period(eq, mime=400.0, eta=0.2,
                                      points_per_period=ppp)
        m = solve_modes_bloch(bg, kx=1.10 * k0, K=0.0, n_modes=1,
                              growth_tol=1e-7)
        return m[0]["gamma"] if m else 0.0

    g32, g64 = gamma_at(32), gamma_at(64)
    assert g32 > 1e-4 and g64 > 1e-4
    assert abs(g32 - g64) / abs(g64) < 1e-2      # 実測 ~1e-14


def test_warm_bloch_eigenfunction_smooth():
    """warm 斜め枝は磁場に結合したモードで、固有関数が滑らか（帯域制限）。"""
    eq = _warm_eq()
    bg = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=48)
    k0 = 2.0 * np.pi / bg.lambda0
    m = solve_modes_bloch(bg, kx=1.10 * k0, K=0.0, n_modes=1, growth_tol=1e-7)[0]
    M = bg.M
    v = m["eigvec"]

    def blk_norm(field: str) -> float:
        i = _field_index(bg.species)[field]
        return float(np.linalg.norm(v[i * M:(i + 1) * M]))

    # δB_z に有意に結合（場と非結合な spurious 表現でない。実測 Bz比 0.065）
    assert blk_norm("Bz") > 0.01 * max(blk_norm(f) for f in _field_index(bg.species))
    # 複素固有ベクトルは fft（rfft は実数専用）。Nyquist 成分が主要成分よりずっと小
    b = _field_index(bg.species)["Bz"]
    bz = v[b * M:(b + 1) * M]
    sp = np.abs(np.fft.fft(bz))
    assert sp[M // 2] < 1e-6 * sp.max()          # Gibbs リンギングなし（実測 3.6e-14）


def _cold_eq() -> dict:
    """cold（PIC対応）SNR 平衡: Te=0.01, Tinc=Tref≈6.25e-4 の鋭いシート。"""
    from fmi.mach_parameters import MachConfig, build_solver_params
    cfg = MachConfig(M_S=154.92, beta_sh=0.25, beta_te=0.1,
                     sigma=0.0, eta=0.2, mime=400.0, gamma=5.0 / 3.0)
    return solve_snr_equilibrium_electron_frame(
        **build_solver_params(cfg), a0_target=0.25, N_points=64, n_steps=25)


def test_fourier_bg_basic():
    """Fourier 背景は n_periods=1・正しい長さ・電荷割当を持ち、密度平均が正規化値。"""
    eq = _warm_eq()
    bg = build_background_fourier(eq, mime=400.0, points_per_period=48)
    assert bg.n_periods == 1 and bg.M == 48
    assert set(bg.species) == {"e", "inc", "ref"}
    assert bg.charge == {"e": -1.0, "inc": 1.0, "ref": 1.0}
    assert bg.lambda0 == pytest.approx(float(eq["lambda0"]))
    assert np.mean(bg.n["e"]) == pytest.approx(1.0, abs=1e-6)          # exp/⟨exp⟩ 正規化
    assert abs(bg.B0z[0]) < 1e-3                                       # 奇対称: y=0 で ~0


def test_fourier_bg_strictly_positive_cold():
    """Fourier 背景は cold の準真空でも密度が厳密に正（barycentric は負に overshoot）。"""
    eq = _cold_eq()
    bg_f = build_background_fourier(eq, mime=400.0, points_per_period=64)
    bg_b = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=64)
    assert min(v.min() for v in bg_f.n.values()) > 0.0                # exp() で正値保証
    assert min(v.min() for v in bg_b.n.values()) < 0.0                # 補間の負値 overshoot


def test_fourier_bg_matches_barycentric_warm_oblique():
    """warm 斜め kx/k0=1.10 の成長率が Fourier 背景と barycentric 背景で一致。"""
    eq = _warm_eq()
    k0 = 2.0 * np.pi / float(eq["lambda0"])
    bg_f = build_background_fourier(eq, mime=400.0, points_per_period=48)
    bg_b = build_background_1period(eq, mime=400.0, eta=0.2, points_per_period=48)
    g_f = solve_modes_bloch(bg_f, kx=1.10 * k0, K=0.0, n_modes=1)[0]["gamma"]
    g_b = solve_modes_bloch(bg_b, kx=1.10 * k0, K=0.0, n_modes=1)[0]["gamma"]
    assert abs(g_f - g_b) / abs(g_b) < 1e-3                           # 実測 ~0%


def test_fourier_bg_cold_dki_converges():
    """cold DKI（kx/k0=0.30, floor=3e-2）成長率が PPP=48→80 で収束（~8.6e-3）。"""
    eq = _cold_eq()
    k0 = 2.0 * np.pi / float(eq["lambda0"])

    def gamma_at(ppp: int) -> float:
        bg = build_background_fourier(eq, mime=400.0, points_per_period=ppp)
        m = solve_modes_bloch(bg, kx=0.30 * k0, K=0.0, n_modes=1,
                              growth_tol=1e-6, density_floor=3e-2)
        return m[0]["gamma"] if m else 0.0

    g48, g80 = gamma_at(48), gamma_at(80)
    assert g48 > 1e-3 and g80 > 1e-3
    assert abs(g48 - g80) / abs(g80) < 0.05                           # 5% 以内で収束


def test_dealias_projector_is_idempotent_and_kills_nyquist():
    """2/3 射影は冪等・実行列で、Nyquist 調波を消し低調波を残す。"""
    from fmi.snr_linear import _dealias_projector

    for M in (48, 96, 300):
        P = _dealias_projector(M)
        assert np.abs(P @ P - P).max() < 1e-12          # 冪等
        assert np.linalg.matrix_rank(P) == 2 * int(2 / 3 * (M / 2)) + 1
        y = np.arange(M) * 2 * np.pi / M
        assert np.abs(P @ np.cos((M // 2) * y)).max() < 1e-12   # Nyquist を消す
        low = np.cos(3 * y)
        assert np.abs(P @ low - low).max() < 1e-12             # 低調波は不変


def test_fourier_diff_matrix_annihilates_nyquist():
    """偶数 M の Fourier 微分行列は Nyquist 調波を消す（偽成長枝の原因）。

    この性質自体は標準だが、そのモードだけが圧力 ∝ k_y^2 と Faraday/Ampère の
    k_y 結合を受けなくなるため、デエイリアシングでの除去が必要になる。
    """
    from fmi.snr_linear import _fourier_diff_matrix

    for M in (48, 96, 300):
        D = _fourier_diff_matrix(M, 2 * np.pi)
        y = np.arange(M) * 2 * np.pi / M
        assert np.abs(D @ np.cos((M // 2) * y)).max() < 1e-10
        assert np.abs(np.linalg.eigvals(D).imag).max() < M // 2   # M/2 に届かない


def test_dealias_preserves_resolved_modes():
    """帯域内に収まっているモードは dealias の有無で変わらない（warm natural）。"""
    from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
    from fmi.snr_linear_bloch import build_background_fourier

    eq = solve_snr_equilibrium_electron_frame(
        eta=0.2, beta_inc=-0.05, beta_ref=0.19, Te=0.5, Tinc=0.1, Tref=0.1,
        a0_target=0.25, N_points=64, n_steps=25)
    k0 = 2.0 * np.pi / float(eq["lambda0"])
    bg = build_background_fourier(eq, mime=400.0, points_per_period=48)
    for r in (0.30, 1.10):
        g_on = np.linalg.eigvals(
            build_operator_bloch(bg, kx=r * k0, K=0.0, dealias=True)).imag.max()
        g_off = np.linalg.eigvals(
            build_operator_bloch(bg, kx=r * k0, K=0.0, dealias=False)).imag.max()
        assert abs(g_on - g_off) / g_off < 1e-6, f"kx/k0={r}: {g_on} vs {g_off}"
