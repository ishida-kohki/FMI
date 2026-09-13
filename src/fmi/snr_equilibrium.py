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


def _clenshaw_curtis_weights(N: int) -> NDArray:
    """Gauss-Lobattoノード上のClenshaw-Curtis積分重みを返す。

    周期にわたる空間平均 <f> = sum(w*f)/sum(w) の計算に使う。
    """
    if N == 0:
        return np.array([1.0])
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    for k in range(N + 1):
        s = 0.0
        for j in range(1, N // 2 + 1):
            b = 2.0 if (2 * j != N) else 1.0
            s += b / (4 * j * j - 1) * np.cos(2 * j * theta[k])
        c = 2.0 if (k != 0 and k != N) else 1.0
        w[k] = c / N * (1.0 - s)
    return w


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
    n_stop: float = 0.0,
    T_stop: float = 0.5,
    u_init: NDArray | None = None,
    _seed_drifts: tuple[float, float, float] | None = None,
) -> dict:
    """3成分SNRプラズマの周期的平衡を解く（任意で静止環境成分を追加）。

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
    n_stop:
        静止環境成分の平均密度（全電子密度 1 に対する比、0 <= n_stop < 1）。
        0（既定）で従来の3成分平衡に厳密に一致する。
        **線形解析の 1/n_{s0} 発散対策としては失敗した**（Notes 参照）。
        既定 0 のまま使わないこと。
    T_stop:
        静止環境成分の温度（電子・イオン共通、m_e*c^2 単位）。
        枯渇の目安は exp(ΔΦ/T_stop) なので、T_stop >~ ΔΦ が必要。
    u_init:
        初期推定ベクトル u=[A(0..N), Phi(0..N), L_half]（長さ 2*N+3）。
        指定時は線形分散 seed と a0 ランプをバイパスし、a0_target で
        1回だけ Newton 解する（フレーム継続の warm-start 用）。
    _seed_drifts:
        線形分散 seed 専用のドリフト (beta_inc, beta_ref, beta_e)。内部用。
        指定時、seed の結合係数 C_AP/C_AA はこの値で計算し、残差・抽出は引数の
        (beta_inc, beta_ref) と そこから導く beta_e を使う。電子静止系（beta_e=0）
        で seed が縮退（Tinc=Tref で C_AP=0）するのを、自然系ドリフトの seed で
        回避して正しい basin を選ぶために solve_snr_equilibrium_electron_frame が使う。

    Returns
    -------
    dict:
        y_full, B0z_full, E0y_full — 1周期分の電磁場プロファイル
        ne_full, ninc_full, nref_full — 密度プロファイル（n̄_e で正規化）
        nstop_e_full, nstop_i_full — 静止環境成分の密度（n_stop=0 なら全て 0）
        lambda0 — 全周期長（c/omega_pe 単位）
        beta_e — 電子ドリフト（電流中性条件から導出）
        params — 入力パラメータの辞書

    Notes
    -----
    **静止環境成分（n_stop > 0）**

    ドリフトしない電子・イオンの対を加える。GEM reconnection challenge
    (Birn et al. 2001, JGR 106, 3715) が Harris 平衡に一様背景 n_b = 0.2 n_0 を
    加えるのと同じ構成で、背景成分は電流を運ばない。

    物理的な位置づけは「周囲プラズマ」というより、単一温度 Maxwellian を仮定した
    流体閉包では表現できない**分布関数の高エネルギー裾**である。実際の分布では
    枯渇領域にも粒子が到達するが、本モデルの n_s ∝ exp(-q_s(Φ-β_s A)/T_s) は
    a0=10 で n_ref を 1e-13〜1e-43 まで落としてしまう。等方な高温サブ成分を
    加えることでこれを補う。

    閉包は他成分と同じだが β=0 なので**ベクトルポテンシャル A が効かない**:

        n_stop_e ∝ exp(+Φ/T_stop),   n_stop_i ∝ exp(-Φ/T_stop)

    既存成分の壊滅的な枯渇は β_s·A の項が作っている（a0=10 で β_ref·ΔA = 2.49 に
    対し ΔΦ = 0.43）ため、静止成分はそれを免れる。T_stop=0.5 なら a0=10 でも
    密度コントラストは 2.4 倍に収まり、実効的な密度下限として機能する。

    規格化は**全電子密度の平均を 1 に保つ**（ω_pe=1 の意味を変えない）:

        n̄_e = 1 - n_stop,  n̄_stop_e = n_stop
        n̄_inc = (1-eta)(1-n_stop),  n̄_ref = eta(1-n_stop),  n̄_stop_i = n_stop

    この配分で平均電荷中性と平均電流中性が同時に満たされ、beta_e の式
    (1-eta)β_inc + eta·β_ref は **n_stop に依存せず不変**である
    （静止成分が電流に寄与せず、(1-n_stop) 因子が両辺で相殺するため）。

    なお本成分は 1/n_ref の発散を解消しない（n_ref 自身は依然として枯渇する）。
    線形解析では snr_linear_relative と併用すること。
    """
    if not 0.0 <= n_stop < 1.0:
        raise ValueError(f"n_stop は 0 <= n_stop < 1 が必要（受領値: {n_stop}）")
    if n_stop > 0.0 and T_stop <= 0.0:
        raise ValueError(f"T_stop は正である必要がある（受領値: {T_stop}）")

    beam = 1.0 - n_stop        # ビーム成分（電子・入射・反射）に配分する密度
    N = N_points
    D_std, x_std = _chebyshev_diff_matrix(N)
    D2_std = D_std @ D_std

    # 周期空間平均 <f> = sum(w*f)/sum(w) 用の積分重み
    cc_w = _clenshaw_curtis_weights(N)
    cc_wsum = cc_w.sum()

    def _avg(f: NDArray) -> float:
        return float(np.sum(cc_w * f) / cc_wsum)

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

        # 密度の前因子は化学ポテンシャル正規化: 空間平均が背景密度に一致するよう
        # exp[...] をその空間平均で割る。これにより周期境界条件が要求する
        # 積分中性条件（平均電荷中性・平均電流中性）が自動的に満たされる。
        e_e = np.exp(arg_e)
        e_i = np.exp(arg_inc)
        e_r = np.exp(arg_ref)
        n_e   = beam         * e_e / _avg(e_e)
        n_inc = beam * (1.0 - eta) * e_i / _avg(e_i)
        n_ref = beam * eta         * e_r / _avg(e_r)

        # 静止環境成分（beta=0 なので A が入らず、電流にも寄与しない）
        if n_stop > 0.0:
            arg_se = np.clip( Phi / T_stop, -100, 100)
            arg_si = np.clip(-Phi / T_stop, -100, 100)
            e_se, e_si = np.exp(arg_se), np.exp(arg_si)
            n_stop_e = n_stop * e_se / _avg(e_se)
            n_stop_i = n_stop * e_si / _avg(e_si)
        else:
            n_stop_e = n_stop_i = 0.0

        # ガウス則（正規化）: d²Phi/dy² = (電子) - (イオン)
        rhs_Phi = (n_e + n_stop_e) - n_inc - n_ref - n_stop_i
        # アンペール則（正規化）: 静止成分は電流を運ばないので寄与しない
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

    if u_init is not None:
        u_init = np.asarray(u_init, dtype=float)
        if u_init.shape != (2 * N + 3,):
            raise ValueError(
                f"u_init の長さは {2*N+3} (=2*N+3) である必要があります。"
                f"実際: {u_init.shape}"
            )
        u_guess = u_init.copy()
        a0_ramp = np.array([a0_target])  # warm-start: 目標振幅で単発solve
    else:
        # --- 線形分散関係に基づく初期推定値 ---
        # seed 専用ドリフト（未指定なら残差と同じ値）。電子静止系で seed 縮退を
        # 避けるため、呼び出し側が自然系ドリフトを渡せるようにする。
        if _seed_drifts is None:
            bi_s, br_s, be_s = beta_inc, beta_ref, beta_e
        else:
            bi_s, br_s, be_s = _seed_drifts
        a0_start = 0.01
        C_PP = 1.0/Te + (1.0-eta)/Tinc + eta/Tref
        C_AP = be_s/Te + (1.0-eta)*bi_s/Tinc + eta*br_s/Tref
        C_AA = -(be_s**2/Te + (1.0-eta)*bi_s**2/Tinc + eta*br_s**2/Tref)

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

    # 密度プロファイル（残差と同じ化学ポテンシャル正規化を適用）
    arg_e   =  (Phi_sol - beta_e   * A_sol) / Te
    arg_inc = -(Phi_sol - beta_inc  * A_sol) / Tinc
    arg_ref = -(Phi_sol - beta_ref  * A_sol) / Tref
    e_e = np.exp(arg_e)
    e_i = np.exp(arg_inc)
    e_r = np.exp(arg_ref)
    n_e_half   = beam         * e_e / _avg(e_e)
    n_inc_half = beam * (1.0 - eta) * e_i / _avg(e_i)
    n_ref_half = beam * eta         * e_r / _avg(e_r)

    if n_stop > 0.0:
        e_se = np.exp(np.clip( Phi_sol / T_stop, -100, 100))
        e_si = np.exp(np.clip(-Phi_sol / T_stop, -100, 100))
        n_stop_e_half = n_stop * e_se / _avg(e_se)
        n_stop_i_half = n_stop * e_si / _avg(e_si)
    else:
        n_stop_e_half = np.zeros_like(n_e_half)
        n_stop_i_half = np.zeros_like(n_e_half)

    # 1周期分に対称拡張:
    #   B, E: 交互フィラメントのため反対称（奇関数）
    #   密度: 対称（偶関数）
    def _mirror_odd(v: NDArray) -> NDArray:
        return np.concatenate([v, -v[-2::-1]])

    def _mirror_even(v: NDArray) -> NDArray:
        return np.concatenate([v, v[-2::-1]])

    return {
        "y_full":    np.concatenate([y_half, 2.0 * L_half - y_half[-2::-1]]),
        "B0z_full":  _mirror_odd(B0z),
        "E0y_full":  _mirror_odd(E0y),
        "ne_full":   _mirror_even(n_e_half),
        "ninc_full": _mirror_even(n_inc_half),
        "nref_full": _mirror_even(n_ref_half),
        "nstop_e_full": _mirror_even(n_stop_e_half),
        "nstop_i_full": _mirror_even(n_stop_i_half),
        "lambda0":   2.0 * L_half,
        "beta_e":    beta_e,
        "u_solution": sol.x.copy(),
        "params": {
            "eta": eta, "beta_inc": beta_inc, "beta_ref": beta_ref,
            "Te": Te, "Tinc": Tinc, "Tref": Tref, "a0_target": a0_target,
            "n_stop": n_stop, "T_stop": T_stop,
        },
    }


def solve_snr_equilibrium_electron_frame(
    eta: float = 0.1,
    beta_inc: float = 0.3,
    beta_ref: float = -0.6,
    Te: float = 0.5,
    Tinc: float = 1.0,
    Tref: float = 0.5,
    a0_target: float = 0.25,
    N_points: int = 64,
    n_steps: int = 25,
    n_stop: float = 0.0,
    T_stop: float = 0.5,
) -> dict:
    """電子静止系（β_e=0）の周期平衡を最初から直接解く。

    支配方程式は solve_snr_equilibrium と同一だが、電子ドリフトを 0 に固定して
    解く。boost や frame-homotopy（自然系で解いてから系を移す2段階）は行わず、
    a0 ランプ1本で電子系平衡を直接求める。

    入力の (beta_inc, beta_ref) は自然系（電流中性系）のドリフトを渡す。内部で
    β_e,nat = (1-eta)*beta_inc + eta*beta_ref を差し引いた電子系ドリフト
    (beta_inc-β_e,nat, beta_ref-β_e,nat) に移す。これらは (1-eta)β_inc'+η β_ref'=0
    を満たし、solve_snr_equilibrium が再導出する beta_e も 0 になる。

    ただし β_e=0 では線形分散 seed の結合係数 C_AP が縮退（Tinc=Tref で C_AP=0）し
    λ0 が崩壊した spurious 解に落ちる。これを避けるため seed のみ自然系ドリフト
    (_seed_drifts) で初期化して正しい basin を選ぶ。残差・抽出（E0y など）は
    β_e=0・電子系ドリフトで行うので、返る E0y は電子系の物理場そのもの。

    Parameters
    ----------
    eta, beta_inc, beta_ref, Te, Tinc, Tref, a0_target, N_points, n_steps,
    n_stop, T_stop:
        solve_snr_equilibrium と同じ。beta_inc/beta_ref は自然系の値。
        静止環境成分は β=0 なので系の取り方に依らず同じ形になる。

    Returns
    -------
    dict:
        solve_snr_equilibrium と同じ形式。beta_e = 0（機械精度）。
        params の beta_inc/beta_ref は電子系ドリフト（自然系値から β_e,nat を減じた値）。
    """
    beta_e_nat = (1.0 - eta) * beta_inc + eta * beta_ref
    bi_e = beta_inc - beta_e_nat
    br_e = beta_ref - beta_e_nat
    return solve_snr_equilibrium(
        eta=eta, beta_inc=bi_e, beta_ref=br_e,
        Te=Te, Tinc=Tinc, Tref=Tref,
        a0_target=a0_target, N_points=N_points, n_steps=n_steps,
        n_stop=n_stop, T_stop=T_stop,
        _seed_drifts=(beta_inc, beta_ref, beta_e_nat),
    )
