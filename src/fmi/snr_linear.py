"""SNR平衡の線形固有値解析（非相対論 warm 多流体 + Maxwell）。

SNR 初期平衡を背景に多流体を線形化し、多周期にタイルした領域で Fourier スペクトル
直接固有値解析を行う。kx をパラメータ、y をスペクトルで解像し、成長率 Im(ω) と
固有ベクトル（モード構造）を求める。導出は docs/snr_linear_modes.md を参照。

正規化は snr_equilibrium と同一（長さ c/ω_pe、β=v/c、c=1、m_e=1, m_i=mime、ω_pe=1）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# 各成分のフィールド名（連続・運動x・運動y）。EM は末尾に Ex, Ey, Bz。
_SPECIES_FIELDS = ("n", "vx", "vy")


@dataclass(frozen=True)
class BackgroundProfiles:
    """線形演算子が要求する多周期背景プロファイル。

    Attributes
    ----------
    y: [0, n_periods*lambda0) 上の一様 M 点。
    n: 成分名 -> n_{s0}(y) 配列（長さ M）。
    beta: 成分名 -> ドリフト u_{s0x}=β_s（定数）。
    T: 成分名 -> 温度 T_s。
    mass: 成分名 -> 質量 m_s。
    charge: 成分名 -> 電荷 q_s。
    B0z: 背景磁場 B0z(y)（長さ M）。
    lambda0: 平衡1周期長。
    n_periods: タイル周期数 N。
    """

    y: NDArray
    n: dict[str, NDArray]
    beta: dict[str, float]
    T: dict[str, float]
    mass: dict[str, float]
    charge: dict[str, float]
    B0z: NDArray
    lambda0: float
    n_periods: int
    E0y: NDArray | None = None

    @property
    def species(self) -> list[str]:
        return list(self.n.keys())

    @property
    def M(self) -> int:
        return self.y.size


def build_background(
    eq: dict,
    mime: float,
    eta: float,
    points_per_period: int = 24,
    n_periods: int = 6,
) -> BackgroundProfiles:
    """平衡ソルバー返り値 eq から多周期背景を作る。

    eq["y_full"]（1周期, 非一様 Chebyshev 鏡映格子）を [0, lambda0) の一様格子へ
    周期補間し、n_periods 回タイルして [0, N*lambda0) を返す。

    Parameters
    ----------
    eq: solve_snr_equilibrium(_electron_frame) の返り値。
    mime: 質量比 m_i/m_e。
    eta: 反射イオン割合（電荷・質量割当には未使用、API 一貫性のため受ける）。
    points_per_period: 1周期あたりの一様格子点数。
    n_periods: タイル周期数。
    """
    lambda0 = float(eq["lambda0"])
    y_src = np.asarray(eq["y_full"], dtype=float)
    p = points_per_period
    y_one = np.linspace(0.0, lambda0, p, endpoint=False)

    def _resample(profile: NDArray) -> NDArray:
        # 周期端を補完してから線形補間（プロファイルは滑らか）
        ys = np.concatenate([y_src, [lambda0]]) #y座標に周期端を追加
        ps = np.concatenate([np.asarray(profile, dtype=float),
                             [np.asarray(profile, dtype=float)[0]]]) #その端点の値を周期性で埋める
        return np.interp(y_one, ys, ps) #線形補間で一様格子へ載せ替える

    ne_one = _resample(eq["ne_full"])
    B0z_one = _resample(eq["B0z_full"])
    n_e = np.tile(ne_one, n_periods)
    n_inc = np.tile(_resample(eq["ninc_full"]), n_periods)
    n_ref = np.tile(_resample(eq["nref_full"]), n_periods)
    B0z = np.tile(B0z_one, n_periods)

    y_full = np.concatenate(
        [y_one + k * lambda0 for k in range(n_periods)]
    )

    params = eq["params"]

    # 背景 E0y は電子の力釣り合い E0y = β_e B0z - T_e n_e'/n_e から格子上で再構成する
    # （Boltzmann 平衡で厳密。鋭い E0y_full を直接補間するより格子と自己無撞着で、
    #  並進ゼロモードを格子誤差の水準で消す）。n_e は滑らかでスペクトル解像が良い。
    D_one = _fourier_diff_matrix(p, lambda0)
    E0y_one = (
        float(eq["beta_e"]) * B0z_one
        - float(params["Te"]) * (D_one @ ne_one) / ne_one
    )
    E0y = np.tile(E0y_one, n_periods)
    return BackgroundProfiles(
        y=y_full,
        n={"e": n_e, "inc": n_inc, "ref": n_ref},
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
        B0z=B0z,
        lambda0=lambda0,
        n_periods=n_periods,
        E0y=E0y,
    )


def _fourier_diff_matrix(M: int, L: float) -> NDArray:
    """周期 [0,L) 一様 M 点上の1階 Fourier 微分行列 D_y（M×M, 実）。

    Trefethen, Spectral Methods in MATLAB の周期微分行列を区間長 L へスケール。
    **偶数 M 専用**（cot 公式）。奇数 M では誤った値を返す
    （M=301 で D @ cos(150y) の振幅が 463、理論値 150 に対し 3 倍ずれる）。

    偶数 M では Nyquist 調波 n=M/2 が消える（D @ cos((M/2)y) = 0）。これは
    Fourier 微分行列の標準的な性質だが、そのモードだけが k_y による安定化を
    受けなくなるため偽の成長枝を生む。`_dealias_projector` で除去すること。
    """
    j = np.arange(M)
    diff = j[:, None] - j[None, :]
    D = np.zeros((M, M))
    mask = diff != 0
    # (pi/L) * (-1)^(j-k) / tan(pi (j-k) / M)
    D[mask] = (np.pi / L) * ((-1.0) ** diff[mask]) / np.tan(np.pi * diff[mask] / M)
    return D


DEALIAS_FRAC = 2.0 / 3.0


def _dealias_projector(M: int, frac: float = DEALIAS_FRAC) -> NDArray:
    """|n| > frac*(M/2) の Fourier 調波を落とす実 M×M 射影行列（2/3 則）。

    2 つの数値障害を同時に除く。

    1. **Nyquist 調波の偽成長**: 偶数 M の Fourier 微分行列は n=M/2 の調波を
       消す（D @ cos((M/2)y) = 0、固有値の虚部も最大 M/2-1 で M/2 に届かない）。
       このモードだけが圧力項 ∝ k_y^2 と Faraday/Ampere の k_y 結合を受けず、
       安定化なしに成長する。SNR 準定常平衡では最速モードのパワーの 97%
       (ppp=300) がこの調波に乗っていた。
    2. **畳み込みのエイリアシング**: 本実装は背景量を対角行列として掛けるため、
       摂動の調波 n と背景の調波 m の積が n+m > M/2 で折り返る。

    P は keep マスクが n→-n で対称なので実行列になる。
        P[j,k] = (1/M) * sum_n keep_n * exp(2*pi*i*n*(j-k)/M)

    Parameters
    ----------
    M: 1周期あたりの格子点数。
    frac: 残す帯域の割合（Nyquist に対する比）。既定 2/3。

    Returns
    -------
    (M, M) の実射影行列（P @ P = P）。
    """
    n = np.fft.fftfreq(M, d=1.0 / M)                  # 整数調波（正負）
    keep = (np.abs(n) <= frac * (M / 2.0)).astype(float)
    phase = np.exp(2j * np.pi * np.outer(np.arange(M), n) / M)
    return np.real((phase * keep) @ phase.conj().T / M)


def _apply_dealias(L: NDArray, P: NDArray, M: int) -> NDArray:
    """ブロック対角射影 kron(I, P) を左右から掛ける（P L P）。"""
    nf = L.shape[0] // M
    A = L.reshape(nf, M, nf, M)
    A = np.einsum("ij,ajbk->aibk", P, A)
    A = np.einsum("ajbk,kl->ajbl", A, P)
    return A.reshape(nf * M, nf * M)


def _field_index(species: list[str]) -> dict[str, int]:
    """フィールド名 -> ブロック番号。各成分 (n,vx,vy) の後に Ex,Ey,Bz。"""
    idx: dict[str, int] = {}
    b = 0
    for s in species:
        for f in _SPECIES_FIELDS:
            idx[f"{s}_{f}"] = b
            b += 1
    for f in ("Ex", "Ey", "Bz"):
        idx[f] = b
        b += 1
    return idx


def build_operator(
    bg: BackgroundProfiles,
    kx: float,
    gamma_ad: float = 5.0 / 3.0,
    coupling: float = 1.0,
    density_floor: float = 0.0,
) -> NDArray:
    """ω x = L x の演算子 L（複素 (3S+3)M × (3S+3)M）を組む。

    各ブロックは docs/snr_linear_modes.md §4 の係数。∂_y は Fourier 微分行列、
    背景 n_{s0}(y),B0z(y) は対角乗算。coupling は電流結合係数（ω_pe=1 で 1）。

    density_floor: 圧力項と背景場力項の 1/n0 を 1/max(n0, density_floor·max n0) に
        置換する相対密度下限。低温平衡で密度が空乏化する領域（n0→0）で 1/n0 が発散し
        非物理な spurious モードを生むのを抑える。0（既定）で従来どおり正則化なし。
        連続・Ampère の n0 乗算（除算でない）は実密度のまま。
        注意: 低温（密度空乏あり）背景では density_floor>0 が実質必須。

    運動 y には背景場が密度揺らぎに及ぼす一次の力 q_s δn_s (E0y - β_s B0z) を含む
    （力釣り合いより T_s δn_s d(ln n_s0)/dy に等しい）。bg.E0y が None の場合は
    E0y=0 として扱う（一様背景の較正テスト用）。
    """
    species = bg.species
    M = bg.M
    S = len(species)
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])  # 周期長 = N*lambda0
    D = _fourier_diff_matrix(M, L_box)
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
        Dn0 = np.diag(n0)
        if density_floor > 0.0:
            n0_reg = np.maximum(n0, density_floor * float(np.max(n0)))
        else:
            n0_reg = n0
        inv_n0 = np.diag(1.0 / n0_reg)
        DB = np.diag(bg.B0z)
        sn, svx, svy = f"{s}_n", f"{s}_vx", f"{s}_vy"

        # 連続
        add(sn, sn, kx * v0 * eye)
        add(sn, svx, kx * Dn0)
        add(sn, svy, -1j * (D @ Dn0))

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
        add(svy, sn, -1j * (gamma_ad * T / m) * (inv_n0 @ D))
        # 背景場が密度揺らぎに及ぼす力 +q_s δn_s (E0y - β_s B0z)
        # （÷(-i m n0) 後は +i(q/m)(E0y - β_s B0z)/n0。並進ゼロモードを消すのに必須）
        add(svy, sn, 1j * (q / m) * np.diag((E0y - v0 * bg.B0z) / n0_reg))

        # Ampère（成分電流の寄与）
        add(Ex, svx, -1j * coupling * q * Dn0)
        add(Ex, sn, -1j * coupling * q * v0 * eye)
        add(Ey, svy, -1j * coupling * q * Dn0)

    # Faraday
    add(Bz, Ey, kx * eye)
    add(Bz, Ex, 1j * D)

    # Ampère（∇×B 項）
    add(Ex, Bz, 1j * D)
    add(Ey, Bz, kx * eye)

    return Lop


def solve_modes(
    bg: BackgroundProfiles,
    kx: float,
    n_modes: int = 6,
    growth_tol: float = 1e-6,
    **op_kw: float,
) -> list[dict]:
    """L を解き、成長モード（Im ω>growth_tol）を成長率降順で返す。

    Returns
    -------
    list of dict: 各要素 {"omega", "gamma", "eigvec"}（成長率降順、最大 n_modes 個）。
    """
    Lop = build_operator(bg, kx, **op_kw)
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


def bloch_wavenumber(mode: dict, bg: BackgroundProfiles) -> float:
    """固有モードの Bloch 波数 K を k0=2π/lambda0 単位で返す（[-0.5, 0.5]）。

    平衡は λ0 周期なので演算子は λ0 並進と可換 → 各固有モードは確定した Bloch K を
    持つ。δB_z を1周期分シフトした並進演算子の固有位相 exp(i K λ0) から K を測る。
    FFT のピーク（格子キャリア）ではなく、この包絡 K が併合（K=k0/2）の指標。
    """
    M = bg.M
    shift = M // bg.n_periods
    bz = mode["eigvec"][-M:]  # 末尾ブロック = Bz
    rho = np.vdot(bz, np.roll(bz, -shift)) / np.vdot(bz, bz)
    return float(np.angle(rho) / (2.0 * np.pi))  # = K / k0


def scan_bloch_spectrum(
    bg: BackgroundProfiles,
    kx: float,
    growth_tol: float = 1e-7,
    **op_kw: float,
) -> dict:
    """成長モードを Bloch セクタ |K|/k0 ごとに集計した分散 γ(K) を返す。

    多周期領域の固有モードは Bloch K = m·k0/N (m=0..N) に整理される。各セクタの
    最大成長率を返す（FMI 分散の Brillouin ゾーン表現）。K=0.5 が併合（period-doubling）。

    Returns
    -------
    dict: {"K_over_k0": NDArray(昇順), "gamma": NDArray} 各セクタの最大 Im(ω)。
    """
    Lop = build_operator(bg, kx, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    N = bg.n_periods
    best: dict[int, float] = {}
    for j in range(eigval.size):
        g = float(eigval[j].imag)
        if g <= growth_tol:
            continue
        K = bloch_wavenumber({"eigvec": eigvec[:, j]}, bg)
        sector = int(round(abs(K) * N))  # 0..N//2
        best[sector] = max(best.get(sector, 0.0), g)
    if not best:
        return {"K_over_k0": np.array([]), "gamma": np.array([])}
    sectors = np.array(sorted(best))
    return {
        "K_over_k0": sectors / N,
        "gamma": np.array([best[s] for s in sectors]),
    }


__all__ = [
    "BackgroundProfiles",
    "build_background",
    "build_operator",
    "solve_modes",
    "bloch_wavenumber",
    "scan_bloch_spectrum",
]
