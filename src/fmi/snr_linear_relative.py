"""相対密度変数 w_s = δn_s / n_{s0} による線形演算子（1/n0 を完全に除去）。

背景
----
``snr_linear_bloch.build_operator_bloch`` は圧力項と背景場力項に 1/n_{s0} を含む。
Harris 的な深い平衡（a0≳2）では反射イオン密度が Boltzmann 的に指数枯渇し、
a0=10 では n_ref,min ≈ 3e-13、すなわち 1/n0 ≈ 3e12 となって演算子の条件数が破綻する。

本モジュールの方針
------------------
δn_s を **相対密度揺らぎ** w_s = δn_s / n_{s0} に取り替える。これは
S = blockdiag(1/n_{s0} on δn blocks, I elsewhere) による**相似変換**
L_rel = S L S^{-1} であり、n_{s0} > 0 である限り**固有値は厳密に不変**である。
密度フロアや環境密度の付加と違い、物理（分散関係）を一切変えない。

変換後、1/n_{s0} は全て ∂_y ln n_{s0} に化ける:

    連続式   (1/n0) ∂_y(n0 δv_y) = ∂_y δv_y + (∂_y ln n0) δv_y
    圧力(y)  (1/n0) ∂_y(δp)      ∝ ∂_y w + (∂_y ln n0) w
    圧力(x)  δn/n0               = w
    背景場力 δn/n0               = w

∂_y ln n_{s0} は指数枯渇プロファイルに対して有界かつ滑らかである
（a0=10 で max|1/n0|=3.5e12 に対し max|∂_y ln n0|=2.17、解像度収束済み）。
逆に Ampère のソース q n0 δv, q v0 n0 w は n0→0 で**小さく**なり、
これは「粒子のいない場所は電流に寄与しない」という正しい極限に対応する。

適用限界
--------
本変換は 1/n_{s0} という**数値**の問題を厳密に取り除くが、密度がほぼゼロの領域で
流体近似が意味を失うという**モデル**の問題は解決しない。むしろ 1/n0 を消したことで、
その枯渇領域の自由度が計算に乗るようになる。

- n_min ≳ 1e-3（warm natural, a0≲2）: 従来演算子・floor と一致。安全。
- n_min ~ 1e-13（eta=0.2, a0=10）: 数値的には健全（max|L|=2.2）だが、成長モード
  上位40個のうち26個が枯渇領域に局在する偽モードになる。
- n_min ~ 1e-43（eta=0.5, T=0.02, a0=10）: 成長モード数が解像度とともに増え
  （ppp 240→400 で 731→1103）、γ_max も 40% 動いて**収束しない**。使用不可。

深く枯渇する平衡を扱うには、枯渇領域を真空にしないこと（静止した準中性の環境
プラズマを平衡の段階で加える）が必要である。詳細は docs/density_regularization.md。

正規化は snr_linear / snr_linear_bloch と同一（c=1, ω_pe=1, m_e=1）。
"""

from __future__ import annotations

import logging

import numpy as np
from numpy.typing import NDArray

from fmi.snr_linear import (
    BackgroundProfiles,
    _apply_dealias,
    _dealias_projector,
    _field_index,
    _fourier_diff_matrix,
)

logger = logging.getLogger(__name__)

__all__ = [
    "build_operator_bloch_relative",
    "solve_modes_bloch_relative",
    "to_absolute_density",
]


def build_operator_bloch_relative(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    gamma_ad: float = 5.0 / 3.0,
    coupling: float = 1.0,
    dealias: float | bool = True,
) -> NDArray:
    """ω x = L_rel(kx,K) x を相対密度変数で組む。

    場の並びは ``_field_index`` と同一だが、``f"{s}_n"`` ブロックの意味が
    δn_s ではなく w_s = δn_s/n_{s0} である点だけが異なる。固有ベクトルを
    従来の δn 表現に戻すには :func:`to_absolute_density` を使う。

    ``build_operator_bloch`` と違い density_floor / n_ambient / viscosity は
    受け取らない。1/n0 が存在しないため正則化パラメータが不要であり、
    結果に恣意的な自由度を残さないためである。

    Args:
        bg: 背景プロファイル（1周期、一様格子）。
        kx: 軸方向波数。
        K: Bloch（横方向）波数。
        dealias: 2/3 デエイリアシング射影を掛けるか。True で 2/3、float で
            その帯域比、False で無効。**既定は True。** 偶数 M の Nyquist 調波は
            Fourier 微分行列に消されて k_y 安定化を受けず偽成長するため
            （`snr_linear._dealias_projector` 参照）、既定で除去する。
        gamma_ad: 断熱指数。
        coupling: Ampère ソースの結合係数（デバッグ用、既定 1）。

    Returns:
        (3S+3)M 次の複素演算子行列。

    Raises:
        ValueError: 背景密度に非正値が含まれる場合（ln n0 が定義できない）。
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

    E0y = bg.E0y if bg.E0y is not None else np.zeros(M)

    for s in species:
        q = bg.charge[s]
        m = bg.mass[s]
        v0 = bg.beta[s]
        T = bg.T[s]
        n0 = bg.n[s]
        if np.any(n0 <= 0.0):
            raise ValueError(
                f"成分 {s} の背景密度に非正値がある（min={float(n0.min()):.3e}）。"
                " 相対密度変数は ln n0 を使うため n0>0 が必要。"
            )
        # 1/n0 の代わりに現れる唯一の背景量。周期関数なので素の D を使う。
        dlogn0 = np.diag(D @ np.log(n0))
        Dn0 = np.diag(n0)
        DB = np.diag(bg.B0z)
        sn, svx, svy = f"{s}_n", f"{s}_vx", f"{s}_vy"

        # 連続（両辺を n0 で割った形）: ω w = kx v0 w + kx δv_x - i(DK + dlnn0)δv_y
        add(sn, sn, kx * v0 * eye)
        add(sn, svx, kx * eye)
        add(sn, svy, -1j * (DK + dlogn0))

        # 運動 x
        add(svx, svx, kx * v0 * eye)
        add(svx, Ex, 1j * q / m * eye)
        add(svx, svy, 1j * q / m * DB)
        add(svx, sn, (kx * gamma_ad * T / m) * eye)

        # 運動 y
        add(svy, svy, kx * v0 * eye)
        add(svy, Ey, 1j * q / m * eye)
        add(svy, svx, -1j * q / m * DB)
        add(svy, Bz, -1j * q * v0 / m * eye)
        add(svy, sn, -1j * (gamma_ad * T / m) * (DK + dlogn0))
        # 背景場が密度揺らぎに及ぼす力 +q_s δn_s (E0y - β_s B0z)/n_{s0}
        add(svy, sn, 1j * (q / m) * np.diag(E0y - v0 * bg.B0z))

        # Ampère（成分電流）: δj = q n0 δv + q v0 δn = q n0 δv + q v0 n0 w
        add(Ex, svx, -1j * coupling * q * Dn0)
        add(Ex, sn, -1j * coupling * q * v0 * Dn0)
        add(Ey, svy, -1j * coupling * q * Dn0)

    # Faraday
    add(Bz, Ey, kx * eye)
    add(Bz, Ex, 1j * DK)

    # Ampère（∇×B）
    add(Ex, Bz, 1j * DK)
    add(Ey, Bz, kx * eye)

    if dealias:
        frac = 2.0 / 3.0 if dealias is True else float(dealias)
        Lop = _apply_dealias(Lop, _dealias_projector(M, frac), M)
    return Lop


def to_absolute_density(vec: NDArray, bg: BackgroundProfiles) -> NDArray:
    """相対密度表現の固有ベクトルを従来の δn 表現へ戻す（w ブロックに n0 を掛ける）。

    パリティ判定や固有関数プロットなど、δn を前提とする既存の解析コードに
    そのまま渡せるようにするための変換。速度・電磁場ブロックは不変。
    """
    out = np.array(vec, dtype=np.complex128, copy=True)
    idx = _field_index(bg.species)
    M = bg.M
    for s in bg.species:
        b = idx[f"{s}_n"]
        out[b * M:(b + 1) * M] *= bg.n[s]
    return out


def solve_modes_bloch_relative(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    n_modes: int = 6,
    growth_tol: float = 1e-6,
    absolute_density: bool = True,
    **op_kw: float,
) -> list[dict]:
    """L_rel(kx,K) を解き、成長モード（Im ω>growth_tol）を成長率降順で返す。

    Args:
        bg: 背景プロファイル。
        kx: 軸方向波数。
        K: Bloch 波数。
        n_modes: 返す最大モード数。
        growth_tol: これ以下の Im ω は成長とみなさない。
        absolute_density: True なら固有ベクトルを δn 表現へ戻して返す
            （既存の解析コードとの互換）。False なら w のまま返す。
        **op_kw: build_operator_bloch_relative へ渡す追加引数。

    Returns:
        {"omega", "gamma", "eigvec"} の辞書のリスト（成長率降順）。
    """
    Lop = build_operator_bloch_relative(bg, kx, K, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    order = np.argsort(-eigval.imag)
    modes: list[dict] = []
    for j in order:
        if eigval[j].imag <= growth_tol:
            break
        v = eigvec[:, j].copy()
        modes.append(
            {
                "omega": complex(eigval[j]),
                "gamma": float(eigval[j].imag),
                "eigvec": to_absolute_density(v, bg) if absolute_density else v,
            }
        )
        if len(modes) >= n_modes:
            break
    return modes
