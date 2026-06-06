"""SNR環境の3成分プラズマ電流フィラメント平衡ソルバー。

Vanthieghem et al. (2018) のペアプラズマ平衡アプローチを
非相対論的な背景電子・入射イオン・反射イオン系に拡張。

正規化: 長さ c/omega_pe、ポテンシャル m_e*c^2/e、温度 m_e*c^2/k_B
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import root


def _chebyshev_diff_matrix(N: int) -> tuple[NDArray, NDArray]:
    """Chebyshev微分行列 D とノード配列 x を返す。"""
    if N == 0:
        return np.array([[0.0]]), np.array([0.0])
    x = np.cos(np.arange(N + 1) * np.pi / N)
    D = np.zeros((N + 1, N + 1))
    p = np.ones(N + 1)
    p[0] = 2.0
    p[N] = 2.0
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = (-1) ** (i + j) * p[i] / (p[j] * (x[i] - x[j]))
    D[0, 0] = (1 + 2 * N**2) / 6.0
    D[N, N] = -(1 + 2 * N**2) / 6.0
    for j in range(1, N):
        D[j, j] = -x[j] / (2 * (1 - x[j] ** 2))
    return D, x


def solve_snr_equilibrium(
    eta: float = 0.1,
    beta_inc: float = 0.3,
    beta_ref: float = -0.6,
    Te: float = 0.5,
    Tinc: float = 1.0,
    Tref: float = 0.5,
    a0_target: float = 0.25,
    N_points: int = 64,
    n_steps: int = 25,
) -> dict:
    """3成分SNRプラズマの周期的平衡を解く。

    Parameters
    ----------
    eta:
        反射イオン割合 n_ref_bar / n_e_bar（0 < eta < 1）。
    beta_inc:
        入射イオンのドリフト速度 / c。
    beta_ref:
        反射イオンのドリフト速度 / c（通常は負）。
    Te:
        電子温度（m_e*c^2 単位）。
    Tinc:
        入射イオン温度（m_e*c^2 単位）。
    Tref:
        反射イオン温度（m_e*c^2 単位）。
    a0_target:
        ベクトルポテンシャル振幅の目標値 e*A_{0x,peak} / (m_e*c^2)。
    N_points:
        Chebyshevの分割数（ノード数は N_points+1）。
    n_steps:
        Newton継続法のステップ数。

    Returns
    -------
    dict:
        y_full, B0z_full, E0y_full — 1周期分の電磁場プロファイル
        ne_full, ninc_full, nref_full — 密度プロファイル（n̄_e で正規化）
        lambda0 — 全周期長（c/omega_pe 単位）
        beta_e — 電子ドリフト（電流中性条件から導出）
        params — 入力パラメータの辞書
    """
    N = N_points
    D_std, x_std = _chebyshev_diff_matrix(N)
    D2_std = D_std @ D_std

    # 電流中性条件から beta_e を導出
    beta_e = (1.0 - eta) * beta_inc + eta * beta_ref

    def _residual(u: NDArray, current_a0: float) -> NDArray:
        A   = u[:N + 1]
        Phi = u[N + 1 : 2 * N + 2]
        L_half = abs(u[-1]) + 1e-6

        dy_dx   = -(L_half / 2.0)
        D_phys  = D_std  / dy_dx
        D2_phys = D2_std / dy_dx**2

        arg_e   = np.clip( (Phi - beta_e   * A) / Te,    -100, 100)
        arg_inc = np.clip(-(Phi - beta_inc  * A) / Tinc,  -100, 100)
        arg_ref = np.clip(-(Phi - beta_ref  * A) / Tref,  -100, 100)

        n_e   = np.exp(arg_e)
        n_inc = (1.0 - eta) * np.exp(arg_inc)
        n_ref = eta          * np.exp(arg_ref)

        # ガウス則（正規化）: d²Phi/dy² = n_e - n_inc - n_ref
        rhs_Phi = n_e - n_inc - n_ref
        # アンペール則（正規化）: d²A/dy² = beta_e*n_e - beta_inc*n_inc - beta_ref*n_ref
        rhs_A   = beta_e * n_e - beta_inc * n_inc - beta_ref * n_ref

        eq_Phi = D2_phys @ Phi - rhs_Phi
        eq_A   = D2_phys @ A   - rhs_A

        factor = dy_dx**2
        res = np.zeros_like(u)
        # 内部ノードの方程式（ノード 1..N-1）
        res[1 : N]           = eq_A[1 : N]   * factor
        res[N + 2 : 2*N + 1] = eq_Phi[1 : N] * factor
        # 境界条件
        res[0]       = (A[0] - current_a0)       * factor   # A(0) = a0
        res[N]       = (D_phys @ A)[N]           * factor   # dA/dy = 0 at y=L_half
        res[N + 1]   = (D_phys @ Phi)[0]         * factor   # dPhi/dy = 0 at y=0
        res[2*N + 1] = (D_phys @ Phi)[N]         * factor   # dPhi/dy = 0 at y=L_half
        res[2*N + 2] = (D_phys @ A)[0]           * factor   # dA/dy = 0 at y=0
        return res

    # --- 線形分散関係に基づく初期推定値 ---
    a0_start = 0.01
    C_PP = 1.0/Te + (1.0-eta)/Tinc + eta/Tref
    C_AP = beta_e/Te + (1.0-eta)*beta_inc/Tinc + eta*beta_ref/Tref
    C_AA = -(beta_e**2/Te + (1.0-eta)*beta_inc**2/Tinc + eta*beta_ref**2/Tref)

    # phi_ratio r が満たす方程式: C_AP*r^2 + (C_AA - C_PP)*r + C_AP = 0
    # （既存コードの a_q=C_APhi, b_q=C_AA-C_PhiPhi, c_q=C_APhi と同じ構造）
    a_q, b_q, c_q = C_AP, C_AA - C_PP, C_AP
    D_q = b_q**2 - 4.0 * a_q * c_q
    if D_q < 0 or abs(C_AP) < 1e-14:
        # 縮退ケース（純静電・ドリフトなし）
        phi_ratio = 0.0
        k2 = -C_AA if -C_AA > 0 else C_PP
    else:
        r1 = (-b_q + np.sqrt(D_q)) / (2.0 * a_q)
        r2 = (-b_q - np.sqrt(D_q)) / (2.0 * a_q)
        phi_ratio, k2 = r1, -(C_AP * r1 + C_AA)
        if k2 <= 0:
            phi_ratio, k2 = r2, -(C_AP * r2 + C_AA)
        if k2 <= 0:
            raise ValueError(
                "指定パラメータで振動解が存在しません。"
                f"C_PP={C_PP:.4f}, C_AP={C_AP:.4f}, C_AA={C_AA:.4f}"
            )

    L_half_guess = np.pi / np.sqrt(k2)
    cos_profile = a0_start * np.cos(np.pi * (1.0 - x_std) / 2.0)
    u_guess = np.zeros(2 * N + 3)
    u_guess[:N + 1]          = cos_profile
    u_guess[N + 1 : 2*N + 2] = phi_ratio * cos_profile
    u_guess[-1]              = L_half_guess

    # --- Newton継続法: a0 を a0_start から a0_target まで段階的に増大 ---
    a0_ramp = np.linspace(a0_start, a0_target, n_steps)
    for i, current_a0 in enumerate(a0_ramp):
        if i > 0:
            u_guess[:2*N + 2] *= current_a0 / a0_ramp[i - 1]
        sol = root(_residual, u_guess, args=(current_a0,), method="lm")
        u_guess = sol.x

    if not sol.success:
        raise RuntimeError(
            f"Newton法がa0={a0_target}で収束しませんでした。"
            f"メッセージ: {sol.message}"
        )

    # --- 解の抽出 ---
    A_sol   = sol.x[:N + 1]
    Phi_sol = sol.x[N + 1 : 2*N + 2]
    L_half  = abs(sol.x[-1])
    D_phys  = D_std / (-(L_half / 2.0))

    B0z = -(D_phys @ A_sol)    # B_{0z} = -dA/dy
    E0y = -(D_phys @ Phi_sol)  # E_{0y} = -dPhi/dy

    # 物理y座標（半周期、c/omega_pe単位）
    y_half = L_half * (1.0 - x_std) / 2.0

    # 密度プロファイル
    arg_e   =  (Phi_sol - beta_e   * A_sol) / Te
    arg_inc = -(Phi_sol - beta_inc  * A_sol) / Tinc
    arg_ref = -(Phi_sol - beta_ref  * A_sol) / Tref
    n_e_half   = np.exp(arg_e)
    n_inc_half = (1.0 - eta) * np.exp(arg_inc)
    n_ref_half = eta          * np.exp(arg_ref)

    # 1周期分に対称拡張:
    #   B, E: 交互フィラメントのため反対称（奇関数）
    #   密度: 対称（偶関数）
    def _mirror_odd(v: NDArray) -> NDArray:
        return np.concatenate([v, -v[-2::-1]])

    def _mirror_even(v: NDArray) -> NDArray:
        return np.concatenate([v, v[-2::-1]])

    return {
        "y_full":    _mirror_even(y_half),
        "B0z_full":  _mirror_odd(B0z),
        "E0y_full":  _mirror_odd(E0y),
        "ne_full":   _mirror_even(n_e_half),
        "ninc_full": _mirror_even(n_inc_half),
        "nref_full": _mirror_even(n_ref_half),
        "lambda0":   2.0 * L_half,
        "beta_e":    beta_e,
        "params": {
            "eta": eta, "beta_inc": beta_inc, "beta_ref": beta_ref,
            "Te": Te, "Tinc": Tinc, "Tref": Tref, "a0_target": a0_target,
        },
    }
