"""SNR平衡の線形固有値解析（Fourier-Bloch 版）。

1周期 [0,λ0) 上で摂動を δ(y)=e^{iKy}û(y)（û は λ0 周期）と分解する。周期背景では
全項の e^{iKy} が消え、演算子は「∂_y → D + iK」の置換だけで得られる。K は明示
パラメータ、kx は軸方向波数。タイル化・Bloch事後抽出・線形補間を用いない。

正規化は snr_linear と同一（c=1, ω_pe=1, m_e=1）。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from fmi.snr_linear import (
    BackgroundProfiles,
    _apply_dealias,
    _dealias_projector,
    _field_index,
    _fourier_diff_matrix,
)


def build_operator_bloch(
    bg: BackgroundProfiles,
    kx: float,
    K: float,
    gamma_ad: float = 5.0 / 3.0,
    coupling: float = 1.0,
    density_floor: float = 0.0,
    n_ambient: float = 0.0,
    viscosity: float = 0.0,
    floor_self_consistent: bool = False,
    dealias: float | bool = True,
) -> NDArray:
    """ω x = L(kx,K) x の演算子を1周期上で組む（∂_y → D + iK）。

    snr_linear.build_operator の写しに Bloch シフト DK=D+iK·I を、∂_y が作用する
    4箇所（連続式 ∂_y(n0δv)、運動yの圧力 ∂_y(δp)、Faraday ∂_yδE_x、Ampère ∇×B）
    にのみ適用する。他ブロックは build_operator と同一。K=0 で build_operator に一致。
    物理式を変更する際は build_operator と本関数の両方を同期させること。

    低密度の谷で発散する 1/n_{s0} の扱いは2通り:

    - ``density_floor`` > 0（build_operator と同一）: 逆数側だけ
      1/max(n_{s0}, floor·max n_{s0}) にクランプする**数値正則化**。順方向の
      Dn0（連続式・Ampère ソース）は素の n_{s0} のままで、自己無撞着ではない。
    - ``n_ambient`` > 0（本関数のみ）: 一様な環境（周囲）プラズマ密度
      n_amb = n_ambient·max_s max_y n_{s0} を全成分に加え、n_{s0}+n_amb を
      **順方向 Dn0 と逆数 1/n_{s0} の両方に一貫適用**する。谷が準真空に達しない
      という物理（SNR フィラメントは周囲プラズマ中にある）に対応し、floor が
      「周囲密度比」という物理パラメータになる。n_ambient>0 のとき density_floor
      は無視する。

      **警告（2026-07-23 検証）**: 現行実装は環境密度を各成分の *ドリフト β_s ごと*
      に加えるため、逆向きドリフトのイオンに一様な対向ビームを注入してしまい、
      **一様な静電二流不安定（γ≈2e-2, kx 非依存, δB_z=0, δn≈1）が混入**する
      （work/diag_ambient_mode.py で確認）。この偽モードが実スペクトルを覆うため、
      n_ambient はこのままでは 1/n0 正則化として使えない。物理的に正しい環境は
      静止（v=0）した準中性の別成分として与える必要がある（未実装）。当面 cold は
      density_floor を使うこと。

    運動 y の背景場力項 q_s δn_s (E0y - β_s B0z) の分母も上と同じ有効密度を使う。

    ``viscosity`` > 0: 各成分の運動方程式(δv_x, δv_y)に人工粘性 +ν·∂_y^2 を追加する
    （Fourier空間で -ν·k_y^2、DK@DK として実装。DKはBlochシフト微分行列なので
    K=0では通常の -ν·k_y^2 拡散に一致）。試験の結果、演算子の非正規性（non-normality）
    により小さなnuでも爆発的な偽成長を誘発することを確認済み（2026-08-06）。
    density_floorの非一貫性を直すものではなく、実用上は推奨しない。デフォルト0で
    既存の全結果に影響しない。

    ``floor_self_consistent`` = True: density_floor>0 のとき、floorされた値
    n0_reg=max(n0,floor·max n0) を **逆数(1/n0)だけでなく順方向Dn0（連続式・
    Ampèreソース）にも同じ値で適用する**（n_ambientと同じ「一貫適用」の考え方だが、
    ドリフトβ_sを持つ環境密度を注入しないぶんシンプル）。谷でのforward/inverseの
    不整合を解消し、高解像度化で悪化するリップル（2026-08-06確認）を抑えることを
    狙う。デフォルトFalseで既存の全結果に影響しない。
    """
    species = bg.species
    M = bg.M
    S = len(species)
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])       # 1周期 = lambda0
    D = _fourier_diff_matrix(M, L_box)
    DK = D + 1j * K * np.eye(M)                   # Bloch シフト
    DKDK = DK @ DK                                 # 人工粘性用 ∂_y^2 相当（-k_y^2）
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

    # 環境密度は全成分共通の絶対量（一様な中性周囲プラズマ）として与える
    n_amb = n_ambient * max(float(np.max(bg.n[s])) for s in species)

    for s in species:
        q = bg.charge[s]
        m = bg.mass[s]
        v0 = bg.beta[s]
        T = bg.T[s]
        n0 = bg.n[s]
        if n_amb > 0.0:
            # 環境密度を順方向・逆数の両方へ一貫適用（自己無撞着な下限）
            n0_eff = n0 + n_amb
            n0_reg = n0_eff
        elif density_floor > 0.0:
            n0_reg = np.maximum(n0, density_floor * float(np.max(n0)))
            # floor_self_consistent: 順方向Dn0にも同じfloor値を使い forward/inverse
            # の不整合を解消する。False（既定）は逆数だけクランプする従来挙動。
            n0_eff = n0_reg if floor_self_consistent else n0
        else:
            n0_eff = n0
            n0_reg = n0
        Dn0 = np.diag(n0_eff)
        inv_n0 = np.diag(1.0 / n0_reg)
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
        if viscosity > 0.0:
            add(svx, svx, 1j * viscosity * DKDK)

        # 運動 y
        add(svy, svy, kx * v0 * eye)
        add(svy, Ey, 1j * q / m * eye)
        add(svy, svx, -1j * q / m * DB)
        add(svy, Bz, -1j * q * v0 / m * eye)
        add(svy, sn, -1j * (gamma_ad * T / m) * (inv_n0 @ DK))
        # 背景場が密度揺らぎに及ぼす力 +q_s δn_s (E0y - β_s B0z)（build_operator と同期）
        add(svy, sn, 1j * (q / m) * np.diag((E0y - v0 * bg.B0z) / n0_reg))
        if viscosity > 0.0:
            add(svy, svy, 1j * viscosity * DKDK)

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

    if dealias:
        frac = 2.0 / 3.0 if dealias is True else float(dealias)
        Lop = _apply_dealias(Lop, _dealias_projector(M, frac), M)
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
    ne_grid = even(ne_h)
    B0z_grid = odd(B_h)
    # 背景 E0y は電子の力釣り合いから格子上で再構成（build_background と同方針）
    D_one = _fourier_diff_matrix(points_per_period, lambda0)
    E0y_grid = (
        float(eq["beta_e"]) * B0z_grid
        - float(params["Te"]) * (D_one @ ne_grid) / ne_grid
    )
    return BackgroundProfiles(
        y=y,
        n={"e": ne_grid, "inc": even(ni_h), "ref": even(nr_h)},
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
        B0z=B0z_grid,
        lambda0=lambda0,
        n_periods=1,
        E0y=E0y_grid,
    )


def _potentials_on_fine_grid(
    eq: dict, M_fine: int
) -> tuple[NDArray, NDArray, NDArray, float]:
    """平衡ポテンシャル A(y),Φ(y) を全周期一様 M_fine 点へスペクトル評価する。

    A,Φ は半周期 Chebyshev 解 eq["u_solution"]=[A(0..N),Φ(0..N),L_half] の先頭
    2(N+1) 成分。両者ともシート中心と周期端で偶対称なので、[0,λ0) の点を [0,L_half]
    へ折り返して Chebyshev barycentric 補間する（build_background_1period と同じ折返し）。
    """
    lam = float(eq["lambda0"])
    L_half = 0.5 * lam
    N = (np.asarray(eq["y_full"]).size - 1) // 2
    u = np.asarray(eq["u_solution"], dtype=float)
    A_half = u[: N + 1]
    Phi_half = u[N + 1 : 2 * N + 2]
    x_cheb, w_cheb = _chebyshev_nodes_weights(N)

    y = np.linspace(0.0, lam, M_fine, endpoint=False)
    s = np.where(y <= L_half, y, lam - y)          # 半周期へ折返し
    xq = 1.0 - 2.0 * s / L_half                    # Chebyshev x へ逆写像
    A_fine = _bary_interp(x_cheb, w_cheb, A_half, xq)     # A は偶
    Phi_fine = _bary_interp(x_cheb, w_cheb, Phi_half, xq)  # Φ は偶
    return y, A_fine, Phi_fine, lam


def build_background_fourier(
    eq: dict,
    mime: float,
    points_per_period: int,
    n_keep: int | None = None,
    oversample: int = 12,
) -> BackgroundProfiles:
    """平衡を Fourier 級数で与える帯域制限背景（正値密度・スペクトル微分）。

    先生方針「Chebyshev 選点を等間隔選点に置換、多めに取って高波数を捨てる」の実装。
    滑らかなポテンシャル A,Φ を細格子でオーバーサンプルして Fourier 変換し、低次
    n_keep 高調波のみ残す（高波数を捨てる）。演算子格子 [0,λ0) の M 点上では:

      - A,Φ とその微分 B0z=-A', E0y=-Φ' を残した Fourier 係数から解析的に評価
        （帯域制限なので格子 D 行列を介さずスペクトル精度で厳密）。
      - 密度 n_{s0}=exp(-q_s(Φ-β_s A)/T_s)/⟨exp⟩ を格子上で再構成 → 正値を保証
        （鋭いシートを直接補間する barycentric 版の負値 overshoot を回避）。

    密度が空乏化する低温平衡では ⟨exp⟩ 正規化後も谷の n0 は極小になる（物理的な準
    真空）。演算子の 1/n0（圧力・背景場力）はこの領域で発散するため、cold では
    build_operator_bloch に density_floor（流体閉包の適用限界）を併用する。

    Parameters
    ----------
    eq: solve_snr_equilibrium(_electron_frame) の返り値。
    mime: 質量比 m_i/m_e。
    points_per_period: 演算子格子の点数 M。
    n_keep: 残す最大 Fourier 高調波次数。None なら M//2-1（格子 Nyquist 未満＝
        エイリアスを出さない最大値）。小さくすると背景をさらに平滑化する。
    oversample: 係数推定用の細格子倍率 M_fine=oversample*M（鋭いシートのエイリアス回避）。
    """
    M = points_per_period
    if n_keep is None:
        n_keep = M // 2 - 1
    n_keep = min(n_keep, M // 2 - 1)               # 格子 Nyquist 未満に制限
    M_fine = oversample * M
    _, A_fine, Phi_fine, lam = _potentials_on_fine_grid(eq, M_fine)

    CA = np.fft.rfft(A_fine)
    CP = np.fft.rfft(Phi_fine)
    K = min(n_keep, CA.size - 1)
    m = np.arange(K + 1)
    km = 2.0 * np.pi * m / lam

    y = np.linspace(0.0, lam, M, endpoint=False)
    phase = np.exp(1j * np.outer(y, km))           # (M, K+1)

    def series(C: NDArray) -> NDArray:
        c = C[: K + 1] / M_fine
        return c[0].real + 2.0 * (phase[:, 1:] @ c[1:]).real

    def dseries(C: NDArray) -> NDArray:
        c = C[: K + 1] / M_fine
        return 2.0 * (phase[:, 1:] @ (1j * km[1:] * c[1:])).real

    A = series(CA)
    Phi = series(CP)
    B0z = -dseries(CA)                             # B0z = -dA/dy（帯域制限で厳密）
    E0y = -dseries(CP)                             # E0y = -dΦ/dy

    p = eq["params"]
    eta = float(p["eta"])
    beta_e = float(eq["beta_e"])
    beta_inc = float(p["beta_inc"])
    beta_ref = float(p["beta_ref"])
    Te, Tinc, Tref = float(p["Te"]), float(p["Tinc"]), float(p["Tref"])

    # 静止環境成分（solve_snr_equilibrium の n_stop>0 のときのみ存在）
    n_stop = float(p.get("n_stop", 0.0))
    T_stop = float(p.get("T_stop", 0.5))
    beam = 1.0 - n_stop

    def density(arg: NDArray, weight: float) -> NDArray:
        e = np.exp(np.clip(arg, -100.0, 100.0))
        return weight * e / e.mean()               # 一様格子なので算術平均=積分平均

    n_e = density((Phi - beta_e * A) / Te, beam)
    n_inc = density(-(Phi - beta_inc * A) / Tinc, beam * (1.0 - eta))
    n_ref = density(-(Phi - beta_ref * A) / Tref, beam * eta)

    n_dict = {"e": n_e, "inc": n_inc, "ref": n_ref}
    beta_dict = {"e": beta_e, "inc": beta_inc, "ref": beta_ref}
    T_dict = {"e": Te, "inc": Tinc, "ref": Tref}
    mass_dict = {"e": 1.0, "inc": float(mime), "ref": float(mime)}
    charge_dict = {"e": -1.0, "inc": 1.0, "ref": 1.0}

    if n_stop > 0.0:
        # β=0 なので A が入らない。既存成分を枯渇させる β_s·A の項を免れる。
        n_dict["stop_e"] = density(Phi / T_stop, n_stop)
        n_dict["stop_i"] = density(-Phi / T_stop, n_stop)
        for name, q, m_s in (("stop_e", -1.0, 1.0), ("stop_i", 1.0, float(mime))):
            beta_dict[name] = 0.0
            T_dict[name] = T_stop
            mass_dict[name] = m_s
            charge_dict[name] = q

    return BackgroundProfiles(
        y=y,
        n=n_dict,
        beta=beta_dict,
        T=T_dict,
        mass=mass_dict,
        charge=charge_dict,
        B0z=B0z,
        lambda0=lam,
        n_periods=1,
        E0y=E0y,
    )


__all__ = [
    "build_operator_bloch",
    "solve_modes_bloch",
    "scan_K_spectrum",
    "reconstruct_eigenfunction",
    "build_background_1period",
    "build_background_fourier",
]
