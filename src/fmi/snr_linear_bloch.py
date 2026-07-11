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
    """ω x = L(kx,K) x の演算子を1周期上で組む（∂_y → D + iK）。"""
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


__all__ = ["build_operator_bloch", "solve_modes_bloch", "scan_K_spectrum"]
