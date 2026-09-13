"""snr_linear.py 解説版（読解専用アノテーション付きコピー）。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
このファイルは src/fmi/snr_linear.py の**内容を等価に保った読解用コピー**です。
実行ロジックは本体と同一で、物理的意図と数値手法を各所に詳しく注釈しています。
本番コードとして import せず（fmi.snr_linear が本線）、コードリーディングの
教材として使ってください。本体を変更したら、この解説版も追随させる想定です。
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【このモジュールが解く問題】
SNR（超新星残骸）前駆体の電流フィラメント平衡を背景として、非相対論の warm 多流体
（電子 e・入射イオン inc・反射イオン ref）に Maxwell 方程式を連立し、線形安定性を
固有値問題として解く。得たいものは:
  - 成長率 γ = Im(ω)  … 各波数でモードが伸びる速さ
  - モード構造（固有ベクトル）… δn, δv, δE, δB の y 依存

【幾何と波数の割り当て】（docs/snr_linear_modes.md と対応）
  - x 軸 = ビーム/電流方向。この向きのゆらぎは平面波 exp(i kx x) で、kx はスカラー
    パラメータとして与える（演算子には kx が数として入る）。
  - y 軸 = フィラメント横断方向。背景 n_{s0}(y), B0z(y) がこの向きに構造を持つので、
    y はスペクトル的に解像する（∂_y を Fourier 微分行列にする）。
  - 背景は周期 λ0（1フィラメント周期）。これを n_periods 回タイルした
    [0, N·λ0) を計算領域にする。周期タイルにより、モードは Bloch 波数 K で分類できる
    （後述の bloch_wavenumber / scan_bloch_spectrum）。

【正規化】（snr_equilibrium と厳密に同一）
  長さ = c/ω_pe、速度 = β = v/c、c = 1、m_e = 1、m_i = mime、ω_pe = 1。
  この規格では固有値 ω は ω_pe 単位、成長率 γ も ω_pe 単位で出る。

【解く固有値問題の形】
  一様背景ドリフト+空間構造をもつ線形化系を、時間 exp(-iωt) で仮定すると
      ω · x = L · x
  という一般化されていない標準固有値問題になる（質量行列は単位）。
  x は全場を y 格子上に並べた長さ (3S+3)·M のベクトル。
  S = 成分数（=3）、M = y 格子点数、各成分が (n, vx, vy) の 3 場、
  さらに電磁場 (Ex, Ey, Bz) の 3 場が末尾に付く。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# 各流体成分がもつ場の名前。順序 (連続=n, 運動x=vx, 運動y=vy) はブロック配置の基準。
# この 3 場 × 成分数 のあとに、共通の電磁場 Ex, Ey, Bz が続く（_field_index 参照）。
_SPECIES_FIELDS = ("n", "vx", "vy")


@dataclass(frozen=True)
class BackgroundProfiles:
    """線形演算子に渡す「凍結した背景平衡」一式（不変データ）。

    frozen=True にしてあるのは、背景は解析中に書き換えない前提を型で保証するため
    （coding-style の Immutability First）。中身の numpy 配列自体は技術的には可変だが、
    このクラスを介して背景を差し替える運用はしない。

    Attributes
    ----------
    y: [0, n_periods*lambda0) 上の一様 M 点。タイル後の全格子。
    n: 成分名 -> n_{s0}(y) 配列（長さ M）。空間構造を持つ唯一の流体量。
    beta: 成分名 -> 背景ドリフト u_{s0x}=β_s（y に依らない定数）。x 方向の流れ。
    T: 成分名 -> 温度 T_s（定数）。圧力 = γ_ad T δn に効く。
    mass: 成分名 -> 質量 m_s（m_e=1, イオン=mime）。
    charge: 成分名 -> 電荷 q_s（電子 -1, イオン +1）。
    B0z: 背景磁場 B0z(y)（長さ M）。フィラメント電流が作る面外磁場。
    lambda0: 平衡1周期長 λ0。
    n_periods: タイル周期数 N。Bloch セクタ数を決める。
    E0y: 背景静電場 E0y(y)。Harris 条件を破る非対称平衡で生じる（後述）。
         None なら 0 扱い（一様背景での較正テスト用）。
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
        # 成分名リスト。dict の挿入順（e, inc, ref）がそのままブロック順になる。
        return list(self.n.keys())

    @property
    def M(self) -> int:
        # y 格子点数。全場ベクトル長 = (3*S + 3) * M の M。
        return self.y.size


def build_background(
    eq: dict,
    mime: float,
    eta: float,
    points_per_period: int = 24,
    n_periods: int = 6,
) -> BackgroundProfiles:
    """平衡ソルバー返り値 eq から、線形解析用の多周期背景を組み立てる。

    やることは2段階:
      (1) 平衡ソルバーが返す非一様 Chebyshev 格子上のプロファイルを、
          [0, λ0) の一様格子へ周期補間で載せ替える（_resample）。
          → Fourier 微分行列は一様格子を要求するため必須の変換。
      (2) 1周期分を n_periods 回タイルして [0, N·λ0) を作る。
          → 周期タイルで Bloch 分類が可能になる。

    Parameters
    ----------
    eq: solve_snr_equilibrium(_electron_frame) の返り値 dict。
        y_full, ne_full, ninc_full, nref_full, B0z_full, lambda0, beta_e, params を含む。
    mime: 質量比 m_i/m_e。イオンの mass に入る。
    eta: 反射イオン割合。ここでは電荷・質量割当に使わず、API 一貫性のため受けるだけ。
    points_per_period: 1周期あたりの一様格子点数 p。分解能パラメータ。
    n_periods: タイル周期数 N。
    """
    lambda0 = float(eq["lambda0"])
    y_src = np.asarray(eq["y_full"], dtype=float)  # 非一様 Chebyshev 鏡映格子（1周期）
    p = points_per_period
    # 目標の一様格子。endpoint=False で [0, λ0) の半開区間（右端 λ0 は含めない=周期点）。
    y_one = np.linspace(0.0, lambda0, p, endpoint=False)

    def _resample(profile: NDArray) -> NDArray:
        """非一様格子上の profile を一様格子 y_one へ周期線形補間する。

        np.interp は範囲外を「端値クランプ」する（外挿しない）ため、そのまま呼ぶと
        y_one の最終点（λ0 近傍）が y_src の右端を超えてクランプされ、周期境界が潰れる。
        そこで x 側に λ0、値側に profile[0]（周期性による端点値）を1点足してから補間する。
        profile が滑らかなので線形補間で十分（warm では誤差小、cold の鋭いシートでは注意）。
        """
        ys = np.concatenate([y_src, [lambda0]])  # y 座標に周期端 λ0 を追加
        ps = np.concatenate([np.asarray(profile, dtype=float),
                             [np.asarray(profile, dtype=float)[0]]])  # 端点値を周期性で埋める
        return np.interp(y_one, ys, ps)  # 一様格子へ線形補間で載せ替え

    # 各流体量を一様格子へ載せ替え、n_periods 回タイル。
    ne_one = _resample(eq["ne_full"])       # 電子密度（1周期, 一様格子）
    B0z_one = _resample(eq["B0z_full"])     # 背景磁場（1周期, 一様格子）
    n_e = np.tile(ne_one, n_periods)                       # 電子密度（全領域）
    n_inc = np.tile(_resample(eq["ninc_full"]), n_periods)  # 入射イオン密度
    n_ref = np.tile(_resample(eq["nref_full"]), n_periods)  # 反射イオン密度
    B0z = np.tile(B0z_one, n_periods)                      # 背景磁場（全領域）

    # 全領域の y 座標。k 番目のタイルは y_one を k·λ0 ずらしたもの。
    y_full = np.concatenate(
        [y_one + k * lambda0 for k in range(n_periods)]
    )

    params = eq["params"]
    
    # 背景静電場 E0y の再構成（ここが設計上のポイント）。
    #   電子の力釣り合い（Boltzmann 平衡）: E0y = β_e B0z - T_e n_e'/n_e。
    # 鋭い E0y_full を直接補間せず、載せ替え後の一様格子上で n_e から再計算する。
    # 理由: 補間した E0y は格子と自己無撞着でなく、背景が線形演算子の「並進ゼロモード」を
    #   厳密に持たなくなる（平衡が方程式を満たさない格子誤差が γ に偽の寄与を出す）。
    #   n_e は滑らかでスペクトル微分の精度が良いので、E0y をこの式で作ると
    #   並進ゼロモードが格子誤差の水準まで消える（数値健全性の担保）。
    D_one = _fourier_diff_matrix(p, lambda0)  # 1周期上の Fourier 微分（n_e' を取るため）
    E0y_one = (
        float(eq["beta_e"]) * B0z_one
        - float(params["Te"]) * (D_one @ ne_one) / ne_one
    )
    E0y = np.tile(E0y_one, n_periods)
    return BackgroundProfiles(
        y=y_full,
        n={"e": n_e, "inc": n_inc, "ref": n_ref},  # dict 順 = ブロック順（e, inc, ref）
        beta={
            "e": float(eq["beta_e"]),          # 電子ドリフト（電子静止系の取り方で決まる）
            "inc": float(params["beta_inc"]),  # 入射イオンドリフト
            "ref": float(params["beta_ref"]),  # 反射イオンドリフト
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
    """周期 [0,L) 一様 M 点上の1階 Fourier 微分行列 D_y（M×M 実行列）を作る。

    「一様周期格子上の値ベクトル」に左から掛けると「その1階 y 微分（格子上）」を返す
    密行列。∂_y をこの行列で表すことで、y 方向を厳密なスペクトル精度で扱える
    （周期・滑らかな背景では有限差分より圧倒的に少ない点数で高精度）。

    出典: Trefethen, "Spectral Methods in MATLAB" の周期微分行列を、標準周期 2π から
    区間長 L へスケール（係数 π/L）。偶数 M を想定した cot 公式だが下式は奇数でも有効。

    行列要素（j≠k）:  D[j,k] = (π/L) · (-1)^(j-k) / tan(π(j-k)/M)
    対角 D[j,j] = 0。
    """
    j = np.arange(M)
    diff = j[:, None] - j[None, :]  # 差 (j-k) の M×M 行列
    D = np.zeros((M, M))
    mask = diff != 0  # 対角以外だけ埋める（対角は 0）
    D[mask] = (np.pi / L) * ((-1.0) ** diff[mask]) / np.tan(np.pi * diff[mask] / M)
    return D


def _field_index(species: list[str]) -> dict[str, int]:
    """フィールド名 -> ブロック番号（0 始まり）の対応表を作る。

    全場ベクトルは M 点ずつのブロックを縦に積んだ構造。ブロック順は:
      [e_n, e_vx, e_vy,  inc_n, inc_vx, inc_vy,  ref_n, ref_vx, ref_vy,  Ex, Ey, Bz]
    各成分の (n, vx, vy) を並べたあと、共通の電磁場 Ex, Ey, Bz を末尾に置く。
    ブロック番号 b の場は全場ベクトルの [b*M : (b+1)*M] を占める。
    """
    idx: dict[str, int] = {}
    b = 0
    for s in species:               # まず各流体成分の 3 場
        for f in _SPECIES_FIELDS:
            idx[f"{s}_{f}"] = b
            b += 1
    for f in ("Ex", "Ey", "Bz"):    # 末尾に電磁場
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
    """固有値問題 ω x = L x の演算子 L を組む（複素 (3S+3)M × (3S+3)M 密行列）。

    各ブロックは docs/snr_linear_modes.md §4 の線形化係数に対応する。実装の対応則:
      - ∂_y      → Fourier 微分行列 D（左から掛ける）
      - 背景量の乗算 n_{s0}(y), B0z(y), 1/n_{s0}(y) → その値を並べた対角行列
      - x 方向の平面波 ∂_x → i·kx（数）。連続・運動の対流項 kx·v0 などに現れる
      - i の因子は exp(-iωt+ikx x) 規約から来る（微分が i を落とす）

    Parameters
    ----------
    kx: x 方向波数（スカラーパラメータ）。
    gamma_ad: 断熱指数 γ_ad（既定 5/3）。圧力項 γ_ad T δn に効く。
    coupling: 電流結合係数。ω_pe=1 規格では 1（Ampère の源電流の重み）。
    density_floor: 圧力項と背景場力項の 1/n0 を 1/max(n0, floor·max n0) に置換する
        相対密度下限。低温平衡で密度が空乏化する領域（n0→0）では 1/n0 が発散し、
        ほぼ真空の格子点が非物理な spurious モードを生む。これを抑える数値正則化。
        0（既定）で正則化なし＝素の 1/n0。
        重要な非対称性: 逆数 1/n0（圧力・背景場力）だけをクランプし、連続・Ampère に
        現れる n0 の「乗算」は実密度のまま。つまり自己無撞着な密度置換ではなく、
        あくまで発散を抑える片側の正則化。
        妥当性の目安（前回の検証結果）:
          - cold の斜め DKI（kx/k0≈0.30）: floor に対して γ が ~5% しか動かない
            → モードは密度が有限の領域に棲む、floor は良性、答えは信頼できる。
          - cold の kx=0 FMI: floor で γ が ~2 倍動く
            → モードが空乏（ほぼ真空）領域に棲み、流体近似が崩れる。運動論が必要。
          - warm（natural）: 密度が空乏化しない（n_min/n_max ≳ 0.12）ので floor は
            そもそも不活性（max(n0, 0.03·max n0)=n0）。natural 解析は floor 非依存。

    背景場の力項について:
      運動 y には、背景場 (E0y, B0z) が密度揺らぎ δn に及ぼす一次の力
      q_s δn_s (E0y - β_s B0z) が入る。力釣り合いより T_s δn_s d(ln n_{s0})/dy に等しく、
      この項が無いと背景の「平行移動」が固有モードにならず、並進ゼロモードが γ に
      偽の寄与を出す。bg.E0y が None のときは E0y=0（一様背景の較正テスト用）。
    """
    species = bg.species
    M = bg.M
    S = len(species)
    # 周期長 = N·λ0。y は半開格子なので、最終点 + 1格子間隔 が全周期長になる。
    L_box = bg.y[-1] + (bg.y[1] - bg.y[0])
    D = _fourier_diff_matrix(M, L_box)  # 全領域上の ∂_y
    eye = np.eye(M)
    idx = _field_index(species)
    dim = (3 * S + 3) * M  # 全場ベクトル長: 各成分 3 場 × S + 電磁 3 場、各 M 点
    Lop = np.zeros((dim, dim), dtype=np.complex128)

    def block(name: str) -> slice:
        # 場名 → 全場ベクトル中のその場の M 点スライス。
        b = idx[name]
        return slice(b * M, (b + 1) * M)

    def add(row: str, col: str, mat: NDArray) -> None:
        # L の (row 場, col 場) ブロックに M×M 行列を加算する糖衣。
        # 「row 場の時間発展方程式に含まれる col 場の係数」を表す。
        Lop[block(row), block(col)] += mat

    Ex, Ey, Bz = "Ex", "Ey", "Bz"

    E0y = bg.E0y if bg.E0y is not None else np.zeros(M)

    # ── 各流体成分ごとに、連続・運動x・運動y・Ampère 源を埋める ──
    for s in species:
        q = bg.charge[s]   # 電荷
        m = bg.mass[s]     # 質量
        v0 = bg.beta[s]    # 背景 x ドリフト β_s
        T = bg.T[s]        # 温度
        n0 = bg.n[s]       # 背景密度 n_{s0}(y)
        Dn0 = np.diag(n0)  # 乗算 ×n0 を表す対角行列
        # 逆数 1/n0 の正則化（density_floor>0 のときだけクランプ）。
        if density_floor > 0.0:
            n0_reg = np.maximum(n0, density_floor * float(np.max(n0)))
        else:
            n0_reg = n0
        inv_n0 = np.diag(1.0 / n0_reg)  # ×(1/n0) を表す対角行列（正則化後）
        DB = np.diag(bg.B0z)            # 乗算 ×B0z(y)
        sn, svx, svy = f"{s}_n", f"{s}_vx", f"{s}_vy"

        # ── 連続の式: ∂_t δn + ∂_x(n0 δvx + v0 δn) + ∂_y(n0 δvy) = 0 ──
        add(sn, sn, kx * v0 * eye)      # 対流 v0 ∂_x δn → kx·v0
        add(sn, svx, kx * Dn0)          # ∂_x(n0 δvx) → kx·n0
        add(sn, svy, -1j * (D @ Dn0))   # ∂_y(n0 δvy) → -i D·n0（i 規約 & 微分）

        # ── 運動 x: ∂_t δvx + v0 ∂_x δvx = (q/m)(δEx + δvy B0z) - (γT/m) ∂_x δn/n0 ──
        add(svx, svx, kx * v0 * eye)             # 対流項
        add(svx, Ex, 1j * q / m * eye)           # 電気力 (q/m) δEx
        add(svx, svy, 1j * q / m * DB)           # ローレンツ力 (q/m) δvy B0z
        add(svx, sn, (kx * gamma_ad * T / m) * inv_n0)  # 圧力勾配 x 成分（1/n0 に floor）

        # ── 運動 y: ∂_t δvy + v0 ∂_x δvy = (q/m)(δEy - δvx B0z - v0 δBz)
        #            - (γT/m) ∂_y δn/n0 + 背景場力 ──
        add(svy, svy, kx * v0 * eye)             # 対流項
        add(svy, Ey, 1j * q / m * eye)           # 電気力 (q/m) δEy
        add(svy, svx, -1j * q / m * DB)          # ローレンツ力 -(q/m) δvx B0z
        add(svy, Bz, -1j * q * v0 / m * eye)     # 背景ドリフト × δBz の力 -(q/m) v0 δBz
        add(svy, sn, -1j * (gamma_ad * T / m) * (inv_n0 @ D))  # 圧力勾配 y 成分
        # 背景場が密度揺らぎに及ぼす力 +q_s δn_s (E0y - β_s B0z)。
        # ÷(-i m n0) 後の形が +i(q/m)(E0y - β_s B0z)/n0。これが並進ゼロモードを消す鍵。
        add(svy, sn, 1j * (q / m) * np.diag((E0y - v0 * bg.B0z) / n0_reg))

        # ── Ampère の源電流 J = Σ q_s (n0 δvs + v0 δn_s) の寄与 ──
        add(Ex, svx, -1j * coupling * q * Dn0)      # x 電流: q n0 δvx
        add(Ex, sn, -1j * coupling * q * v0 * eye)  # x 電流: q v0 δn（ドリフト対流分）
        add(Ey, svy, -1j * coupling * q * Dn0)      # y 電流: q n0 δvy

    # ── Faraday: ∂_t δBz = -(∂_x δEy - ∂_y δEx) ──
    add(Bz, Ey, kx * eye)   # -∂_x δEy → kx（符号は i 規約込み）
    add(Bz, Ex, 1j * D)     # +∂_y δEx → i D

    # ── Ampère の ∇×B 項: ∂_t δE = ... + c²(∇×B)（c=1）──
    add(Ex, Bz, 1j * D)     # (∇×B)_x = ∂_y Bz → i D
    add(Ey, Bz, kx * eye)   # (∇×B)_y = -∂_x Bz → kx

    return Lop


def solve_modes(
    bg: BackgroundProfiles,
    kx: float,
    n_modes: int = 6,
    growth_tol: float = 1e-6,
    **op_kw: float,
) -> list[dict]:
    """L を作って対角化し、成長モード（Im ω > growth_tol）を成長率降順で返す。

    密行列の完全対角化（np.linalg.eig）で全固有値を求め、成長するものだけ拾う。
    **op_kw は build_operator にそのまま渡る（gamma_ad, coupling, density_floor など）。

    Returns
    -------
    list of dict: 各要素 {"omega": 複素固有値, "gamma": Im(ω), "eigvec": 固有ベクトル}。
        成長率降順、最大 n_modes 個。成長モードが無ければ空リスト。
    """
    Lop = build_operator(bg, kx, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    order = np.argsort(-eigval.imag)  # Im(ω) 降順（成長率が大きい順）
    modes: list[dict] = []
    for j in order:
        if eigval[j].imag <= growth_tol:
            break  # 降順なので、閾値を切ったら以降は全て非成長 → 打ち切り
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
    """固有モードの Bloch 波数 K を k0=2π/λ0 単位で返す（値域 [-0.5, 0.5]）。

    背景は λ0 周期なので、線形演算子 L は「λ0 の空間並進」と可換。よって各固有モードは
    確定した Bloch 位相 exp(i K λ0) を持つ（Bloch/Floquet の定理）。この K を、δBz を
    1周期分ずらす並進演算子の期待位相から測る。

    実装: δBz を shift=M/N（=1周期の格子点数）だけ巡回シフトし、元との重なり
      ρ = <bz | roll(bz, -shift)> / <bz|bz> ≈ exp(i K λ0) を作る。
      その偏角を 2π で割ると K/k0 が出る。
    注意: FFT のピーク（格子キャリア波数）ではなく、この「包絡の Bloch K」が
      period-doubling / フィラメント併合（K=k0/2）の指標になる。
    """
    M = bg.M
    shift = M // bg.n_periods  # 1周期分の格子点数
    bz = mode["eigvec"][-M:]   # 全場ベクトル末尾ブロック = δBz
    rho = np.vdot(bz, np.roll(bz, -shift)) / np.vdot(bz, bz)  # ≈ exp(i K λ0)
    return float(np.angle(rho) / (2.0 * np.pi))  # = K / k0


def scan_bloch_spectrum(
    bg: BackgroundProfiles,
    kx: float,
    growth_tol: float = 1e-7,
    **op_kw: float,
) -> dict:
    """成長モードを Bloch セクタ |K|/k0 ごとに集計した分散 γ(K) を返す。

    多周期領域（N タイル）の固有モードは、Bloch K = m·k0/N（m=0..N）の離散セクタに
    整理される。各セクタで最大の成長率を拾って返す＝FMI 分散の Brillouin ゾーン表現。
    K=0（全周期同位相）は基本フィラメントモード、K=0.5 は併合（period-doubling）に対応。

    Returns
    -------
    dict: {"K_over_k0": NDArray（昇順の各セクタ値）, "gamma": NDArray（各セクタ最大 Im ω）}。
        成長モードが無ければ両方空配列。
    """
    Lop = build_operator(bg, kx, **op_kw)
    eigval, eigvec = np.linalg.eig(Lop)
    N = bg.n_periods
    best: dict[int, float] = {}  # セクタ番号 -> そのセクタの最大成長率
    for j in range(eigval.size):
        g = float(eigval[j].imag)
        if g <= growth_tol:
            continue
        K = bloch_wavenumber({"eigvec": eigvec[:, j]}, bg)
        sector = int(round(abs(K) * N))  # |K|/k0·N を丸めてセクタ 0..N//2 に量子化
        best[sector] = max(best.get(sector, 0.0), g)
    if not best:
        return {"K_over_k0": np.array([]), "gamma": np.array([])}
    sectors = np.array(sorted(best))
    return {
        "K_over_k0": sectors / N,  # セクタ番号を k0 単位の Bloch 波数へ戻す
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
