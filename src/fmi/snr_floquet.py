"""SNR 3成分平衡の Floquet（モノドロミー）線形解析。

既存のタイル法 Fourier ソルバ（snr_linear）とは独立の並行トラック。
Vanthieghem+2018 の Floquet 手法（1周期のモノドロミー積分）を SNR 3成分・
非相対論モデルに適用する。平衡は Chebyshev スペクトル解（snr_equilibrium）を
barycentric 補間で任意 y で厳密評価する（一様格子への再サンプルをしない）。

Step 1（本ファイルの現段階）: 平衡→背景関数 bg_func(y) の橋渡し。
  bg_func(y) は任意 y で n_s0(y), B0z(y), E0y(y), d ln n_s0/dy を返す。
  勾配は力釣り合い d ln n_{s0}/dy = q_s (E0y - v_{s0x} B0z)/T_s から代数的に
  与える（数値微分なし）。以降の N 行列（∂_y x = N x）はこれらの値のみで組める。
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .snr_equilibrium import _chebyshev_diff_matrix, _clenshaw_curtis_weights

_SPECIES = ("e", "inc", "ref")


@dataclass(frozen=True)
class BackgroundPoint:
    """bg_func(y) の返り値: ある y での背景量（3成分）。"""

    n: NDArray        # n_{s0}(y)          shape (3,)
    B0z: float        # B_{0z}(y)
    E0y: float        # E_{0y}(y)
    dln_n: NDArray    # d ln n_{s0}/dy     shape (3,)  （力釣り合いから）
    v0: NDArray       # v_{s0x}=β_s (一定) shape (3,)
    q: NDArray        # 電荷 (-1,+1,+1)    shape (3,)
    T: NDArray        # 温度               shape (3,)


def _bary_weights(x: NDArray) -> NDArray:
    """標準 Chebyshev-Gauss-Lobatto ノード用の barycentric 重み w_j=(-1)^j, 端点半分。"""
    n = x.size - 1
    w = np.ones(n + 1)
    w[1::2] = -1.0
    w[0] *= 0.5
    w[-1] *= 0.5
    return w


def _bary_eval(xn: NDArray, w: NDArray, fn: NDArray, xq: float) -> float:
    """barycentric 補間。xn 上の値 fn を点 xq で評価（xq がノード上なら厳密値）。"""
    d = xq - xn
    hit = np.isclose(d, 0.0)
    if np.any(hit):
        return float(fn[np.argmax(hit)])
    c = w / d
    return float(c @ fn / c.sum())


def build_bg_func(eq: dict) -> tuple[Callable[[float], BackgroundPoint], float]:
    """平衡 eq から bg_func(y) と周期 lambda0 を作る。

    eq は solve_snr_equilibrium(_electron_frame) の返り値。u_solution=[A(0..N),
    Phi(0..N), L_half] を生の Chebyshev 解として使い、半周期 [0,L_half] の
    Chebyshev ノード上で A,Phi,B0z,E0y を評価、barycentric 補間で任意 y を張る。
    y は鏡映で [0,lambda0=2 L_half) をカバー（A,Phi,n=偶、B0z,E0y=奇）。
    """
    u = np.asarray(eq["u_solution"], dtype=float)
    N = (u.size - 3) // 2
    A = u[: N + 1]
    Phi = u[N + 1 : 2 * N + 2]
    L_half = abs(u[-1])
    lambda0 = 2.0 * L_half

    D_std, x_std = _chebyshev_diff_matrix(N)
    D_phys = D_std / (-(L_half / 2.0))
    B0z_nodes = -(D_phys @ A)      # B_{0z} = -dA/dy（半周期ノード上）
    E0y_nodes = -(D_phys @ Phi)    # E_{0y} = -dPhi/dy
    w = _bary_weights(x_std)

    p = eq["params"]
    beta = np.array([float(eq["beta_e"]), float(p["beta_inc"]), float(p["beta_ref"])])
    T = np.array([float(p["Te"]), float(p["Tinc"]), float(p["Tref"])])
    q = np.array([-1.0, 1.0, 1.0])
    eta = float(p["eta"])
    frac = np.array([1.0, 1.0 - eta, eta])  # n_s0 の種別スケール

    # 密度正規化 norm_s = <exp(arg_s)>（平衡と同じ化学ポテンシャル正規化）。
    # arg_s = -q_s (Phi - beta_s A)/T_s をノード上で作り、_avg と同じ台形平均で規格化。
    cc_w = _clenshaw_curtis_weights(N)
    cc_wsum = cc_w.sum()

    def _avg(v: NDArray) -> float:
        # snr_equilibrium._avg と同一: Clenshaw-Curtis 加重平均
        return float(np.sum(cc_w * v) / cc_wsum)

    norm = np.empty(3)
    for i in range(3):
        arg = -q[i] * (Phi - beta[i] * A) / T[i]
        norm[i] = _avg(np.exp(arg))

    def bg_func(y: float) -> BackgroundPoint:
        yy = float(y) % lambda0
        if yy <= L_half:
            yf, sgn = yy, 1.0
        else:
            yf, sgn = lambda0 - yy, -1.0   # 鏡映（B0z,E0y は奇→sgn=-1）
        xq = 1.0 - 2.0 * yf / L_half
        A_y = _bary_eval(x_std, w, A, xq)
        Phi_y = _bary_eval(x_std, w, Phi, xq)
        B0z = sgn * _bary_eval(x_std, w, B0z_nodes, xq)
        E0y = sgn * _bary_eval(x_std, w, E0y_nodes, xq)
        n = np.empty(3)
        for i in range(3):
            arg = -q[i] * (Phi_y - beta[i] * A_y) / T[i]
            n[i] = frac[i] * np.exp(arg) / norm[i]
        # 力釣り合い: d ln n_{s0}/dy = q_s (E0y - v_{s0x} B0z)/T_s
        dln = q * (E0y - beta * B0z) / T
        return BackgroundPoint(n=n, B0z=B0z, E0y=E0y, dln_n=dln,
                               v0=beta, q=q, T=T)

    return bg_func, lambda0


# ---------------------------------------------------------------------------
# Step 2: N 行列 (∂_y x = N(y,ω,kx) x, 8次元)
# 状態 x = [δE_x, δB_z, δn_e, δn_inc, δn_ref, δv_ey, δv_incy, δv_refy]
# δv_sx（運動x, 代数）と δE_y（Ampère-y, 代数）を消去。導出は docs 参照。
# ---------------------------------------------------------------------------

_DIM = 8


def get_N_matrix_snr(
    y: float,
    omega: complex,
    kx: float,
    bg_func: Callable[[float], BackgroundPoint],
    gamma_ad: float,
    mass: NDArray,
) -> NDArray:
    """1階系 ∂_y x = N x の 8×8 複素行列 N(y,ω,kx)。"""
    bp = bg_func(y)
    n, B0z, E0y = bp.n, bp.B0z, bp.E0y
    v0, q, T, dln = bp.v0, bp.q, bp.T, bp.dln_n
    m = mass
    Om = omega - kx * v0            # Ω_s = ω - kx β_s  (shape 3)

    N = np.zeros((_DIM, _DIM), dtype=np.complex128)
    iE, iB = 0, 1

    def i_n(s: int) -> int:
        return 2 + s

    def i_v(s: int) -> int:
        return 5 + s

    # --- ∂_y δE_x  (Faraday + δE_y 消去) ---
    N[iE, iB] = 1j * (kx**2 / omega - omega)
    for s in range(3):
        N[iE, i_v(s)] = (kx / omega) * q[s] * n[s]

    # --- ∂_y δB_z  (Ampère-x + δv_sx 消去) ---
    sum_Ex = np.sum(q**2 * n / (m * Om))
    N[iB, iE] = 1j * (sum_Ex - omega)
    for s in range(3):
        N[iB, i_n(s)] = q[s] * kx * gamma_ad * T[s] / (m[s] * Om[s]) + q[s] * v0[s]
        N[iB, i_v(s)] = 1j * q[s]**2 * n[s] * B0z / (m[s] * Om[s])

    # --- ∂_y δn_s  (運動y + δv_sx, δE_y 消去)。両辺を γ_ad T_s で割る ---
    for s in range(3):
        r = i_n(s)
        gT = gamma_ad * T[s]
        N[r, iE] = (-1j * q[s]**2 * n[s] * B0z / (m[s] * Om[s])) / gT
        N[r, iB] = (q[s] * n[s] * (kx / omega - v0[s])) / gT
        N[r, i_n(s)] += (-q[s] * B0z * kx * gamma_ad * T[s] / (m[s] * Om[s])
                         + q[s] * (E0y - v0[s] * B0z)) / gT
        for sp in range(3):
            coeff = -1j * q[s] * n[s] * q[sp] * n[sp] / omega       # δE_y 経由（全 s'）
            if sp == s:
                coeff += (-1j * q[s]**2 * n[s] * B0z**2 / (m[s] * Om[s])
                          + 1j * m[s] * n[s] * Om[s])
            N[r, i_v(sp)] += coeff / gT

    # --- ∂_y δv_sy  (連続 + δv_sx 消去) ---
    for s in range(3):
        r = i_v(s)
        N[r, iE] = kx * q[s] / (m[s] * Om[s])
        N[r, i_n(s)] = (1j * omega / n[s] - 1j * kx * v0[s] / n[s]
                        - 1j * kx**2 * gamma_ad * T[s] / (m[s] * n[s] * Om[s]))
        N[r, i_v(s)] = kx * q[s] * B0z / (m[s] * Om[s]) - dln[s]

    return N


def monodromy_multipliers(
    omega: complex,
    kx: float,
    bg_func: Callable[[float], BackgroundPoint],
    lambda0: float,
    mass: NDArray,
    gamma_ad: float = 5.0 / 3.0,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> NDArray:
    """1周期 [0,λ0] を積分して基本行列 M(λ0)（=モノドロミー, M(0)=I）の
    固有値＝特性乗数 σ を返す。Bloch 波数は K=arg(σ)/λ0。"""
    from scipy.integrate import solve_ivp

    def rhs(y: float, Mf: NDArray) -> NDArray:
        M = (Mf[:_DIM**2] + 1j * Mf[_DIM**2:]).reshape(_DIM, _DIM)
        Nm = get_N_matrix_snr(y, omega, kx, bg_func, gamma_ad, mass)
        dM = Nm @ M
        return np.concatenate([dM.real.ravel(), dM.imag.ravel()])

    M0 = np.eye(_DIM, dtype=np.complex128)
    M0f = np.concatenate([M0.real.ravel(), M0.imag.ravel()])
    sol = solve_ivp(rhs, [0.0, lambda0], M0f, method="Radau", rtol=rtol, atol=atol)
    Mend = (sol.y[:_DIM**2, -1] + 1j * sol.y[_DIM**2:, -1]).reshape(_DIM, _DIM)
    return np.linalg.eigvals(Mend)
