"""SNR平衡の線形固有値解析（Fourier-Bloch 版）。

1周期 [0,λ0) 上で摂動を δ(y)=e^{iKy}û(y)（û は λ0 周期）と分解する。周期背景では
全項の e^{iKy} が消え、演算子は「∂_y → D + iK」の置換だけで得られる。K は明示
パラメータ、kx は軸方向波数。タイル化・Bloch事後抽出・線形補間を用いない。

正規化は snr_linear と同一（c=1, ω_pe=1, m_e=1）。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from fmi.snr_linear import BackgroundProfiles, _field_index, _fourier_diff_matrix


def build_operator_bloch(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    gamma_ad: float = 5.0 / 3.0,
    coupling: float = 1.0,
) -> NDArray:
    """ω x = L(kx,K) x の演算子を1周期上で組む（∂_y → D + iK）。

    snr_linear.build_operator の写しに Bloch シフト DK=D+iK·I を、∂_y が作用する
    4箇所（連続式 ∂_y(n0δv)、運動yの圧力 ∂_y(δp)、Faraday ∂_yδE_x、Ampère ∇×B）
    にのみ適用する。他ブロックは build_operator と同一。K=0 で build_operator に一致。
    物理式を変更する際は build_operator と本関数の両方を同期させること。
    """
    species = bg.species
    M = bg.M
    S = len(species)
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])       # 1周期 = lambda0
    D = _fourier_diff_matrix(M, L_box)
    DK = D + 1j * K * np.eye(M)                   # Bloch シフト
    eye = np.eye(M)
    idx = _field_index(species)
    dim = (3 * S + 3) * M
    Lop = np.zeros((dim, dim), dtype=np.complex128)

    def block(name: str) -> slice:
        b = idx[name]
        return slice(b * M, (b + 1) * M)

    def add(row: str, col: str, mat: NDArray) -> None:
        Lop[block(row), block(col)] += mat

    Ex, Ey, Bz = "Ex", "Ey", "Bz"

    for s in species:
        q = bg.charge[s]
        m = bg.mass[s]
        v0 = bg.beta[s]
        T = bg.T[s]
        n0 = bg.n[s]
        Dn0 = np.diag(n0)
        inv_n0 = np.diag(1.0 / n0)
        DB = np.diag(bg.B0z)
        sn, svx, svy = f"{s}_n", f"{s}_vx", f"{s}_vy"

        # 連続
        add(sn, sn, kx * v0 * eye)
        add(sn, svx, kx * Dn0)
        add(sn, svy, -1j * (DK @ Dn0))

        # 運動 x
        add(svx, svx, kx * v0 * eye)
        add(svx, Ex, 1j * q / m * eye)
        add(svx, svy, 1j * q / m * DB)
        add(svx, sn, (kx * gamma_ad * T / m) * inv_n0)

        # 運動 y
        add(svy, svy, kx * v0 * eye)
        add(svy, Ey, 1j * q / m * eye)
        add(svy, svx, -1j * q / m * DB)
        add(svy, Bz, -1j * q * v0 / m * eye)
        add(svy, sn, -1j * (gamma_ad * T / m) * (inv_n0 @ DK))

        # Ampère（成分電流）
        add(Ex, svx, -1j * coupling * q * Dn0)
        add(Ex, sn, -1j * coupling * q * v0 * eye)
        add(Ey, svy, -1j * coupling * q * Dn0)

    # Faraday
    add(Bz, Ey, kx * eye)
    add(Bz, Ex, 1j * DK)

    # Ampère（∇×B）
    add(Ex, Bz, 1j * DK)
    add(Ey, Bz, kx * eye)

    return Lop


def solve_modes_bloch(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    n_modes: int = 6,
    growth_tol: float = 1e-6,
    **op_kw: float,
) -> list[dict]:
    """L(kx,K) を解き、成長モード（Im ω>growth_tol）を成長率降順で返す。"""
    Lop = build_operator_bloch(bg, kx, K, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    order = np.argsort(-eigval.imag)
    modes: list[dict] = []
    for j in order:
        if eigval[j].imag <= growth_tol:
            break
        modes.append(
            {
                "omega": complex(eigval[j]),
                "gamma": float(eigval[j].imag),
                "eigvec": eigvec[:, j].copy(),
            }
        )
        if len(modes) >= n_modes:
            break
    return modes


def scan_K_spectrum(
    bg: BackgroundProfiles,
    kx: float,
    K_over_k0,
    growth_tol: float = 1e-7,
    **op_kw: float,
) -> dict:
    """Bloch 波数 K を明示走査し、各 K の最大成長率 γ(K) を返す。"""
    k0 = 2.0 * np.pi / bg.lambda0
    r = np.asarray(K_over_k0, dtype=float)
    gamma = np.zeros(r.size)
    for i, ri in enumerate(r):
        modes = solve_modes_bloch(bg, kx, ri * k0, n_modes=1,
                                  growth_tol=growth_tol, **op_kw)
        gamma[i] = modes[0]["gamma"] if modes else 0.0
    return {"K_over_k0": r, "gamma": gamma}


def _chebyshev_nodes_weights(N: int) -> tuple[NDArray, NDArray]:
    """Chebyshev-Gauss-Lobatto ノード x=cos(jπ/N)（降順 1→-1）と barycentric 重み。"""
    x = np.cos(np.arange(N + 1) * np.pi / N)
    w = np.ones(N + 1)
    w[1::2] = -1.0
    w[0] *= 0.5
    w[N] *= 0.5
    return x, w


def _bary_interp(x_nodes: NDArray, w: NDArray, f_nodes: NDArray,
                 x_query) -> NDArray:
    """Berrut–Trefethen barycentric 補間（ノード一致時は値を直接返す）。"""
    xq = np.atleast_1d(np.asarray(x_query, dtype=float))
    out = np.empty(xq.shape, dtype=float)
    for i, xx in enumerate(xq):
        diff = xx - x_nodes
        hit = np.isclose(diff, 0.0)
        if np.any(hit):
            out[i] = f_nodes[int(np.argmax(hit))]
        else:
            t = w / diff
            out[i] = float((t @ f_nodes) / t.sum())
    return out


def reconstruct_eigenfunction(
    mode: dict,
    bg: BackgroundProfiles,
    K: float,
    field: str = "Bz",
    n_display: int = 6,
) -> tuple[NDArray, NDArray]:
    """固有ベクトル û から δ(y)=e^{iKy}û(y) を n_display 周期展開して返す。

    Args:
        mode: solve_modes_bloch が返すモード辞書（キー "eigvec"）
        bg: 背景プロファイル
        K: Bloch 波数
        field: フィールド名（デフォルト "Bz"）
        n_display: 表示用に何周期分展開するか（デフォルト 6）

    Returns:
        (ys, delta): 拡張 y 座標と δ(y)=e^{iKy}û(y)
    """
    M = bg.M
    lam = bg.lambda0
    b = _field_index(bg.species)[field]
    uhat = mode["eigvec"][b * M:(b + 1) * M]
    ys = np.concatenate([bg.y + k * lam for k in range(n_display)])
    delta = np.exp(1j * K * ys) * np.tile(uhat, n_display)
    return ys, delta


def build_background_1period(
    eq: dict,
    mime: float,
    eta: float,
    points_per_period: int,
) -> BackgroundProfiles:
    """平衡（Chebyshev 半周期解）を1周期一様格子へ barycentric スペクトル補間で載せる。

    平衡の *_full は半周期 Chebyshev 解の鏡映（B0z は奇、密度は偶）。先頭 N+1 点が
    半周期値。これを Chebyshev barycentric で任意 y に評価し、鏡映パリティを適用して
    [0,λ0) の一様格子を作る（np.interp の C⁰ 折れ点を排除）。

    `eta` は密度分割済みの `eq` には未使用（tiled 版 build_background との API 一貫性
    のために受ける）。

    注意: スペクトル補間は正値性を保存しない（ノード間で節点最小値を下回る over/under-
    shoot が起こりうる）。warm 平衡（密度最小 ~0.1、帯域制限）では無害だが、密度が 0 に
    迫る鋭い平衡では運動yの 1/n0 が発散しうる（cold σ=0 の既知限界。plan 参照）。
    """
    lambda0 = float(eq["lambda0"])
    L_half = 0.5 * lambda0
    N = (np.asarray(eq["y_full"]).size - 1) // 2
    x_cheb, w_cheb = _chebyshev_nodes_weights(N)

    def half(name: str) -> NDArray:
        return np.asarray(eq[name], dtype=float)[: N + 1]

    ne_h, ni_h, nr_h = half("ne_full"), half("ninc_full"), half("nref_full")
    B_h = half("B0z_full")

    y = np.linspace(0.0, lambda0, points_per_period, endpoint=False)
    s = np.where(y <= L_half, y, lambda0 - y)     # 半周期へ折り返し
    xq = 1.0 - 2.0 * s / L_half                    # y_half=L_half(1-x)/2 の逆写像

    def even(vals_h: NDArray) -> NDArray:
        return _bary_interp(x_cheb, w_cheb, vals_h, xq)

    def odd(vals_h: NDArray) -> NDArray:
        sign = np.where(y <= L_half, 1.0, -1.0)
        return sign * _bary_interp(x_cheb, w_cheb, vals_h, xq)

    params = eq["params"]
    return BackgroundProfiles(
        y=y,
        n={"e": even(ne_h), "inc": even(ni_h), "ref": even(nr_h)},
        beta={
            "e": float(eq["beta_e"]),
            "inc": float(params["beta_inc"]),
            "ref": float(params["beta_ref"]),
        },
        T={
            "e": float(params["Te"]),
            "inc": float(params["Tinc"]),
            "ref": float(params["Tref"]),
        },
        mass={"e": 1.0, "inc": float(mime), "ref": float(mime)},
        charge={"e": -1.0, "inc": 1.0, "ref": 1.0},
        B0z=odd(B_h),
        lambda0=lambda0,
        n_periods=1,
    )


__all__ = [
    "build_operator_bloch",
    "solve_modes_bloch",
    "scan_K_spectrum",
    "reconstruct_eigenfunction",
    "build_background_1period",
]
