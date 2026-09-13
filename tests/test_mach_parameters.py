"""マッハ数基準パラメータ設計器の回帰テスト。

基準は現行 PIC ラン config_2d_inp_xy.toml:
    sigma=0.0025, mime=400, ush=0.125, alpha=0.2, betae=8, betai=betar=0.5
これが M_A=100, M_S≈155 に対応し、Te=0.01, Tinc=6.25e-4 を再現することを確認する。
"""

import pytest

from fmi.mach_parameters import (
    MachConfig,
    build_pic_config,
    build_solver_params,
    to_electron_rest_frame,
)

CURRENT = MachConfig(
    M_A=100.0,
    M_S=155.0,
    eta=0.2,
    mime=400.0,
    sigma=0.0025,
    gamma=5.0 / 3.0,
    beta_e_plasma=8.0,
)


def test_solver_params_reproduce_current_run():
    sp = build_solver_params(CURRENT)
    # M_S=155 は丸め値（betai=0.5 の厳密値は M_S=154.92）なので Tinc は 0.1% ずれる。
    assert sp["Tinc"] == pytest.approx(6.25e-4, rel=2e-3)
    assert sp["Te"] == pytest.approx(0.01, rel=1e-6)  # Te は M_S 非依存（厳密）
    assert sp["Tref"] == pytest.approx(sp["Tinc"])


def test_exact_M_S_gives_betai_half():
    # betai=0.5 を厳密に与える M_S を使えば Tinc は厳密に 6.25e-4
    import math
    M_S_exact = CURRENT.M_A / math.sqrt(CURRENT.gamma * 0.5 / 2.0)
    sp = build_solver_params(MachConfig(
        M_A=CURRENT.M_A, M_S=M_S_exact, eta=CURRENT.eta, mime=CURRENT.mime,
        sigma=CURRENT.sigma, gamma=CURRENT.gamma, beta_e_plasma=8.0,
    ))
    assert sp["Tinc"] == pytest.approx(6.25e-4, rel=1e-9)


def test_current_neutrality_beta_e_near_zero():
    sp = build_solver_params(CURRENT)
    beta_e = (1 - CURRENT.eta) * sp["beta_inc"] + CURRENT.eta * sp["beta_ref"]
    assert abs(beta_e) < 1e-2


def test_pic_config_roundtrip_betas():
    pic = build_pic_config(CURRENT)
    assert pic["betae"] == pytest.approx(8.0, rel=1e-2)
    assert pic["betai"] == pytest.approx(0.5, rel=1e-2)
    assert pic["ush"] == pytest.approx(0.125, rel=2e-2)
    assert pic["needs_direct_temperature"] is False


def test_beta_plasma_matches_mach_identity():
    # β_i = (2/γ)(M_A/M_S)²
    pic = build_pic_config(CURRENT)
    expect = (2.0 / CURRENT.gamma) * (CURRENT.M_A / CURRENT.M_S) ** 2
    assert pic["betai"] == pytest.approx(expect, rel=1e-6)


def test_unmagnetized_requires_velocity_anchor():
    bad = MachConfig(M_A=float("nan"), M_S=20.0, sigma=0.0, beta_te=0.1)
    with pytest.raises(ValueError):
        build_solver_params(bad)  # beta_sh 未指定でエラー
    good = MachConfig(
        M_A=float("nan"), M_S=20.0, sigma=0.0, beta_sh=0.25, beta_te=0.1
    )
    sp = build_solver_params(good)
    pic = build_pic_config(good)
    assert sp["Tinc"] > 0
    assert pic["needs_direct_temperature"] is True


def test_unmagnetized_form_reproduces_magnetized_params():
    """本プロジェクト既定の非磁化(σ=0)指定が磁化版と同一の β・T を出すこと。

    σ は平衡 ODE に現れないため、速度 Vsh と電子温度 vte を直接与えれば
    磁化版と機械精度で一致する（HTML/docstring の主張のロック）。
    """
    import math
    sigma = CURRENT.sigma
    beta_sh = CURRENT.M_A * math.sqrt(sigma / CURRENT.mime)  # = 0.25
    beta_te = math.sqrt(sigma * CURRENT.beta_e_plasma / 2.0)  # = 0.1
    unmag = MachConfig(
        M_S=CURRENT.M_S, beta_sh=beta_sh, beta_te=beta_te,
        sigma=0.0, eta=CURRENT.eta, mime=CURRENT.mime, gamma=CURRENT.gamma,
    )
    a = build_solver_params(CURRENT)
    b = build_solver_params(unmag)
    for key in a:
        assert b[key] == pytest.approx(a[key], rel=1e-12, abs=1e-15)


def test_electron_rest_frame_zeroes_beta_e():
    sp = build_solver_params(CURRENT)
    erf = to_electron_rest_frame(sp)
    # 電子静止系では電流中性条件の電子ドリフトが厳密に 0
    beta_e = (1 - CURRENT.eta) * erf["beta_inc"] + CURRENT.eta * erf["beta_ref"]
    assert abs(beta_e) < 1e-9


def test_electron_rest_frame_preserves_temperatures_and_eta():
    sp = build_solver_params(CURRENT)
    erf = to_electron_rest_frame(sp)
    # 温度と eta は不変、ドリフトのみブースト
    assert erf["eta"] == sp["eta"]
    assert erf["Te"] == sp["Te"]
    assert erf["Tinc"] == sp["Tinc"]
    assert erf["Tref"] == sp["Tref"]
    # ブースト量 ≈ 元の β_e（小さい）
    beta_e0 = (1 - CURRENT.eta) * sp["beta_inc"] + CURRENT.eta * sp["beta_ref"]
    assert erf["beta_inc"] - sp["beta_inc"] == pytest.approx(-beta_e0, rel=1e-2)


def test_electron_rest_frame_relativistic_branch_residual_is_order_beta2():
    sp = build_solver_params(CURRENT)
    erf = to_electron_rest_frame(sp, relativistic=True)
    # 相対論ブーストでは solver が再導出する β_e に O(β²) の残差が残る:
    # 厳密 0 ではないが小さく、非相対論ブーストより緩い。
    beta_e = (1 - CURRENT.eta) * erf["beta_inc"] + CURRENT.eta * erf["beta_ref"]
    assert 0.0 < abs(beta_e) < 1e-4


def test_electron_knob_exclusive():
    # 2つ同時指定はエラー
    two = MachConfig(M_A=100, M_S=155, sigma=0.0025, beta_e_plasma=8.0, beta_te=0.1)
    with pytest.raises(ValueError):
        build_solver_params(two)
    # 0個指定もエラー
    none = MachConfig(M_A=100, M_S=155, sigma=0.0025)
    with pytest.raises(ValueError):
        build_solver_params(none)
