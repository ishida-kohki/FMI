import numpy as np
import pytest
from fmi.snr_equilibrium import (
    solve_snr_equilibrium,
    solve_snr_equilibrium_electron_frame,
)


def test_import():
    assert callable(solve_snr_equilibrium)


def test_snr_default_converges():
    """デフォルトSNRパラメータでエラーなく収束すること。"""
    result = solve_snr_equilibrium(
        eta=0.1, beta_inc=0.3, beta_ref=-0.6,
        Te=0.5, Tinc=1.0, Tref=0.5,
        a0_target=0.1,   # 高速テストのため小振幅
        N_points=32, n_steps=10,
    )
    assert "lambda0" in result
    assert result["lambda0"] > 0


def test_current_neutrality_derived():
    """beta_e == (1-eta)*beta_inc + eta*beta_ref が成立すること。"""
    eta, beta_inc, beta_ref = 0.1, 0.3, -0.6
    result = solve_snr_equilibrium(
        eta=eta, beta_inc=beta_inc, beta_ref=beta_ref,
        Te=0.5, Tinc=1.0, Tref=0.5,
        a0_target=0.1, N_points=32, n_steps=10,
    )
    beta_e_expected = (1 - eta) * beta_inc + eta * beta_ref
    assert abs(result["beta_e"] - beta_e_expected) < 1e-12


def test_B_antisymmetric():
    """B0z_full は反対称であること（B[0] == -B[-1]）。"""
    result = solve_snr_equilibrium(
        eta=0.1, beta_inc=0.3, beta_ref=-0.6,
        Te=0.5, Tinc=1.0, Tref=0.5,
        a0_target=0.1, N_points=32, n_steps=10,
    )
    B = result["B0z_full"]
    assert abs(B[0] + B[-1]) < 1e-6, f"B の反対称性破れ: B[0]={B[0]:.6f}, B[-1]={B[-1]:.6f}"


def test_field_BC_satisfied_at_fold():
    """折り返し点（y=0, y=L_half）で磁場・電場がゼロになること。

    密度の化学ポテンシャル正規化により積分中性条件が満たされ、
    Neumann境界条件 dA/dy=dPhi/dy=0 が機械精度で成立する。
    これを破ると折り返し点に不連続（キンク）が現れる。
    """
    N = 64
    result = solve_snr_equilibrium(
        eta=0.1, beta_inc=0.3, beta_ref=-0.6,
        Te=0.5, Tinc=1.0, Tref=0.5,
        a0_target=0.25, N_points=N, n_steps=25,
    )
    B = result["B0z_full"]
    E = result["E0y_full"]
    # 半周期端点 = 全体配列の index 0（y=0）と index N（y=L_half）
    assert abs(B[0]) < 1e-8, f"B0z(y=0) が非ゼロ: {B[0]:.2e}"
    assert abs(B[N]) < 1e-8, f"B0z(y=L_half) が非ゼロ: {B[N]:.2e}"
    assert abs(E[0]) < 1e-8, f"E0y(y=0) が非ゼロ: {E[0]:.2e}"
    assert abs(E[N]) < 1e-8, f"E0y(y=L_half) が非ゼロ: {E[N]:.2e}"
