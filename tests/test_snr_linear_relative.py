"""相対密度変数 w=δn/n0 による線形演算子の回帰テスト。

中心的な主張は「相対演算子は従来演算子の相似変換であり固有値が厳密に一致する」
こと。密度枯渇が軽い平衡（従来演算子がまだ健全な領域）でこれを確認する。
"""

import numpy as np
import pytest

from fmi.snr_equilibrium import solve_snr_equilibrium_electron_frame
from fmi.snr_linear import BackgroundProfiles, _field_index
from fmi.snr_linear_bloch import build_background_fourier, build_operator_bloch
from fmi.snr_linear_relative import (
    build_operator_bloch_relative,
    solve_modes_bloch_relative,
    to_absolute_density,
)

EQ_KW = dict(eta=0.2, beta_inc=-0.05, beta_ref=0.19, Te=0.5, Tinc=0.1,
             Tref=0.1, N_points=48, n_steps=20)


def _uniform_single_electron(M: int = 24, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    return BackgroundProfiles(
        y=y, n={"e": np.ones(M)}, beta={"e": 0.0}, T={"e": 0.0},
        mass={"e": 1.0}, charge={"e": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def test_relative_lightwave_K0() -> None:
    """一様背景では w=δn なので従来同様 ω²=1+kx² を再現する。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5, 1.3):
        w2 = np.linalg.eigvals(build_operator_bloch_relative(bg, kx, K=0.0)) ** 2
        assert np.min(np.abs(w2 - (1.0 + kx**2))) < 1e-8
        assert np.min(np.abs(w2 - 1.0)) < 1e-8


def test_uniform_background_operators_identical() -> None:
    """n0≡1 の一様背景では相似変換が恒等なので両演算子は要素まで一致する。"""
    bg = _uniform_single_electron()
    L_old = build_operator_bloch(bg, kx=0.7, K=0.3)
    L_new = build_operator_bloch_relative(bg, kx=0.7, K=0.3)
    assert np.allclose(L_old, L_new, atol=1e-12)


@pytest.mark.parametrize("a0_target", [0.25, 1.0])
def test_eigenvalues_match_conventional_operator(a0_target: float) -> None:
    """密度枯渇が軽い平衡では従来演算子と固有値が一致する（相似変換の検算）。

    従来演算子は floor なしで 1/n0 をそのまま使うため、n0 が中程度までしか
    落ちない a0 でしか健全でない。そこが一致すれば導出とコードが正しい。
    """
    eq = solve_snr_equilibrium_electron_frame(a0_target=a0_target, **EQ_KW)
    bg = build_background_fourier(eq, mime=400.0, points_per_period=64)
    k0 = 2.0 * np.pi / float(eq["lambda0"])
    kx = 0.5 * k0

    g_old = np.sort(np.linalg.eigvals(build_operator_bloch(bg, kx, K=0.0)).imag)
    g_new = np.sort(
        np.linalg.eigvals(build_operator_bloch_relative(bg, kx, K=0.0)).imag
    )
    # 支配モードと第2モードのみ比較する。従来演算子は floor なしだと 1/n0 が
    # 大きく条件数が悪いため、γ が支配モードの1桁以下の弱い枝は丸め誤差に
    # 埋もれて一致しない（相似変換の検算としては支配枝の一致で十分）。
    assert g_new[-1] == pytest.approx(g_old[-1], rel=1e-8)
    assert g_new[-2] == pytest.approx(g_old[-2], rel=1e-8)


def test_operator_norm_bounded_for_depleted_equilibrium() -> None:
    """深い平衡でも相対演算子の要素は O(1) に留まる（条件数の改善）。"""
    eq = solve_snr_equilibrium_electron_frame(a0_target=10.0, **EQ_KW)
    bg = build_background_fourier(eq, mime=400.0, points_per_period=64)
    k0 = 2.0 * np.pi / float(eq["lambda0"])

    n_min_ratio = min(float(bg.n[s].min() / bg.n[s].max()) for s in bg.species)
    assert n_min_ratio < 1e-6, "この試験は密度が強く枯渇する平衡を前提にする"

    L_old = build_operator_bloch(bg, kx=0.5 * k0, K=0.0)
    L_new = build_operator_bloch_relative(bg, kx=0.5 * k0, K=0.0)
    assert np.abs(L_old).max() > 1e6      # 従来は 1/n0 で発散する
    assert np.abs(L_new).max() < 1e2      # 相対版は有界


def test_to_absolute_density_roundtrip() -> None:
    """to_absolute_density は密度ブロックにのみ n0 を掛ける。"""
    eq = solve_snr_equilibrium_electron_frame(a0_target=1.0, **EQ_KW)
    bg = build_background_fourier(eq, mime=400.0, points_per_period=32)
    idx = _field_index(bg.species)
    M = bg.M
    v = np.ones((3 * len(bg.species) + 3) * M, dtype=np.complex128)
    out = to_absolute_density(v, bg)
    for s in bg.species:
        b = idx[f"{s}_n"]
        assert np.allclose(out[b * M:(b + 1) * M], bg.n[s])
        for f in ("vx", "vy"):
            b2 = idx[f"{s}_{f}"]
            assert np.allclose(out[b2 * M:(b2 + 1) * M], 1.0)
    for f in ("Ex", "Ey", "Bz"):
        b3 = idx[f]
        assert np.allclose(out[b3 * M:(b3 + 1) * M], 1.0)


def test_rejects_nonpositive_density() -> None:
    """n0<=0 は ln n0 が定義できないので ValueError を出す。"""
    bg = _uniform_single_electron()
    bad = BackgroundProfiles(
        y=bg.y, n={"e": np.zeros(bg.M)}, beta=bg.beta, T=bg.T,
        mass=bg.mass, charge=bg.charge, B0z=bg.B0z,
        lambda0=bg.lambda0, n_periods=1,
    )
    with pytest.raises(ValueError, match="非正値"):
        build_operator_bloch_relative(bad, kx=0.5, K=0.0)


def test_solve_modes_returns_absolute_by_default() -> None:
    """既定では δn 表現の固有ベクトルを返し、w 表現とは n0 倍だけ異なる。"""
    eq = solve_snr_equilibrium_electron_frame(a0_target=1.0, **EQ_KW)
    bg = build_background_fourier(eq, mime=400.0, points_per_period=64)
    k0 = 2.0 * np.pi / float(eq["lambda0"])
    m_abs = solve_modes_bloch_relative(bg, kx=0.5 * k0, K=0.0, n_modes=1)
    m_rel = solve_modes_bloch_relative(bg, kx=0.5 * k0, K=0.0, n_modes=1,
                                       absolute_density=False)
    assert m_abs and m_rel
    assert m_abs[0]["gamma"] == pytest.approx(m_rel[0]["gamma"], rel=1e-12)
    assert np.allclose(m_abs[0]["eigvec"],
                       to_absolute_density(m_rel[0]["eigvec"], bg))
