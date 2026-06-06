import numpy as np
import pytest
from fmi.snr_equilibrium import solve_snr_equilibrium


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
