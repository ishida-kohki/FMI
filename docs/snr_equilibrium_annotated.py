"""snr_equilibrium.py 解説版（読解専用アノテーション付きコピー）。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
このファイルは src/fmi/snr_equilibrium.py の**内容を等価に保った読解用コピー**です。
実行ロジックは本体と同一で、物理的意図と数値手法を各所に詳しく注釈しています。
本番コードとして import せず（fmi.snr_equilibrium が本線）、コードリーディングの
教材として使ってください。本体を変更したら、この解説版も追随させる想定です。
（対応する図解ドキュメント: docs/snr_equilibrium_solver.html）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【このモジュールが解く問題】
SNR（超新星残骸）衝撃波前面近傍に生じる、電子(e)・入射イオン(inc)・反射イオン(ref)
の3成分プラズマが作る**電流フィラメントの周期的非線形平衡解**を求める。
Vanthieghem et al. (2018, Phys. Plasmas 25, 072115) の相対論的ペアプラズマ平衡
（Chebyshevスペクトル法+Newton継続法）を、非相対論3成分系へ拡張したもの。

【3成分プラズマ系】（η∈(0,1) は反射イオン割合）
  成分        電荷    平均密度目安        ドリフト    温度
  背景電子 e   -1      n̄_e（正規化基準=1）  β_e         T_e
  入射イオン   +1      (1-η) n̄_e           β_inc       T_inc
  反射イオン   +1      η n̄_e               β_ref       T_ref

【正規化】（snr_linear と厳密に同一）
  長さ = c/ω_pe、速度 = β = v/c（c=1）、スカラーポテンシャル Φ = eφ/(m_e c²)、
  ベクトルポテンシャル A = eA_x/(m_e c²)、温度 T_s = k_B T_s^phys/(m_e c²)。
  質量比 m_i/m_e はこの平衡方程式には陽に現れない
  （snr_linear の線形安定性解析で初めて必要になる）。

【密度分布（Boltzmann平衡）】
  各成分が局所熱平衡にあると仮定し、非相対論極限の1粒子エネルギー
  （静電ポテンシャルエネルギー q_sΦ と、ドリフト運動が canonical momentum を
  通じて感じるベクトルポテンシャルの効果 -q_s β_s A）から、密度は
      n_s/n̄_s ∝ exp[ -q_s(Φ - β_s A)/T_s ]
  比例係数（化学ポテンシャル正規化定数）は、指数因子の**空間平均**が
  背景密度に一致するように決める（詳細は _residual 内のコメント参照）。
  電子 q_e=-1 なので arg_e = +(Φ-β_e A)/Te、イオン q=+1 なので
  arg_inc/ref = -(Φ-β_{inc,ref} A)/T_{inc,ref} と符号が反転する。

【支配方程式】（1次元、y方向のみに空間構造）
  Poisson（ガウス則）:  d²Φ/dy² = n_e - n_inc - n_ref
  Ampere:              d²A/dy² = β_e n_e - β_inc n_inc - β_ref n_ref
  電荷密度 ρ = -n_e+n_inc+n_ref、電流密度 J_x = -β_e n_e+β_inc n_inc+β_ref n_ref
  とすると、上式はそれぞれ d²Φ/dy² = -ρ、d²A/dy² = -J_x という
  Heaviside-Lorentz単位系の標準形（B0z=-dA/dy, E0y=-dΦ/dy の符号規約と整合）。

【中性条件】
  電荷中性: n̄_inc=(1-η)n̄_e, n̄_ref=ηn̄_e の定義から自動的に成立（追加式不要）。
  電流中性: x方向全電流ゼロ ⟹ β_e=(1-η)β_inc+ηβ_ref。これは入力パラメータから
  解析的に求まるので、Newton の未知数には含めない（1行の代入で済む）。

【境界条件と半周期ドメイン】
  解は周期λ0=2·L_half の中で対称性を持つので、半周期 ỹ∈[0, L_half] だけを解く。
  L_half（半周期長）はNewton系の未知数の一部（固有値的に決まる）。
    y=0       : A(0)=a0（振幅固定）, dA/dy=0（B0z=0, 反射イオンフィラメント中心の対称性）,
                dΦ/dy=0（E0y=0, 同じ対称性）
    y=L_half  : dA/dy=0（B0z=0, 入射イオンフィラメント中心=逆向き電流の対称面）,
                dΦ/dy=0（E0y=0, 同じ対称面）
  5個の境界条件 + 2(N-1)個の内部方程式 = 2N+3個 = 未知数の数（A:N+1, Φ:N+1, L_half:1）
  としてちょうど閉じた正方系になる。

【数値手法】
  1. Chebyshev-Gauss-Lobattoノード上でスペクトル微分行列を構成（_chebyshev_diff_matrix）。
  2. 密度の空間平均には Clenshaw-Curtis 求積重みを使う（_clenshaw_curtis_weights）。
  3. 小振幅極限の線形分散関係から初期推定値 (A,Φ,L_half) を解析的に構成。
  4. Newton継続法（ホモトピー）: a0 を a0_start=0.01 から a0_target まで
     段階的に増やしながら scipy.optimize.root（method="lm"）で逐次収束させる。
  5. 半周期解を対称拡張（B0z,E0y は奇関数、密度は偶関数）して1周期分の解を返す。

【ペアプラズマとの違い（Vanthieghem 2018 との対比）】
  2成分ペアプラズマでは電子・陽電子の指数関数の差が sinh にまとまり、
  支配方程式が sinh 型 Liouville 方程式に帰着し、積分中性条件は対称性から
  自動的に満たされる。3成分系ではこの対称性がないため、指数因子を空間平均で
  正規化する処理が本質的に必要（_residual 内の n_e, n_inc, n_ref の式を参照）。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import root


def _chebyshev_diff_matrix(N: int) -> tuple[NDArray, NDArray]:
    """Chebyshev微分行列 D とノード配列 x を返す（標準区間 x∈[-1,1]）。

    Chebyshev-Gauss-Lobattoノード x_k = cos(kπ/N), k=0,...,N（N+1点、両端を含む）
    上で、関数値ベクトル f = [f(x_0),...,f(x_N)] に D を左から掛けると
    その1階微分 f'(x_k) の近似が返る（密行列、スペクトル精度）。
    等間隔格子の有限差分と違い、滑らかな関数に対して指数関数的に速く収束する
    （これが平衡ソルバーの主要な離散化手法）。

    行列要素（Trefethen, "Spectral Methods in MATLAB" の公式そのもの）:
      D[i,j] = (-1)^(i+j) p_i/(p_j (x_i-x_j))   (i≠j; p_0=p_N=2, それ以外1)
      D[0,0] = (1+2N²)/6、D[N,N] = -(1+2N²)/6   （両端の特別な対角成分）
      D[j,j] = -x_j / (2(1-x_j²))               （内部ノードの対角成分）

    N=0 の退化ケースは 1×1 のゼロ行列（微分不要な自明系）を返す。
    """
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

    周期にわたる空間平均 <f> = sum(w*f)/sum(w) の計算に使う（_avg 関数参照）。
    Chebyshev微分と同じノード上で定義された求積法なので、平衡ソルバーが使う
    離散化（微分・平均とも同一グリッド）が内部で一致し、数値誤差が
    機械精度レベルに抑えられる（詳細は _residual 内コメント参照）。

    導出: Chebyshev係数のコサイン級数を項別積分した閉形式（標準的な公式）。
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
    u_init: NDArray | None = None,
    _seed_drifts: tuple[float, float, float] | None = None,
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
        lambda0 — 全周期長（c/omega_pe 単位）
        beta_e — 電子ドリフト（電流中性条件から導出）
        params — 入力パラメータの辞書
    """
    N = N_points
    # 標準区間 x∈[-1,1] 上の Chebyshev 微分行列。D2_std は2階微分（D_std を2回掛けた行列）。
    D_std, x_std = _chebyshev_diff_matrix(N)
    D2_std = D_std @ D_std

    # 周期空間平均 <f> = sum(w*f)/sum(w) 用の積分重み
    cc_w = _clenshaw_curtis_weights(N)
    cc_wsum = cc_w.sum()

    def _avg(f: NDArray) -> float:
        return float(np.sum(cc_w * f) / cc_wsum)

    # 電流中性条件から beta_e を導出（x方向の全電流がゼロになるよう電子ドリフトを固定）。
    # これは Newton の未知数ではなく、入力パラメータから解析的に決まる定数。
    beta_e = (1.0 - eta) * beta_inc + eta * beta_ref

    def _residual(u: NDArray, current_a0: float) -> NDArray:
        """Newton法に渡す残差ベクトル R(u)。R(u)=0 の解が平衡解。

        u のレイアウト: [A(0),...,A(N), Phi(0),...,Phi(N), L_half]（長さ 2N+3）。
        res も同じ長さで、内部ノードの物理方程式（Poisson/Ampere）と
        境界条件（Neumann×4、振幅固定×1）を1本のベクトルに詰める。
        """
        A   = u[:N + 1]
        Phi = u[N + 1 : 2 * N + 2]
        # L_half（半周期の物理長）は符号自由なので絶対値を取り、0除算回避に微小値を足す。
        L_half = abs(u[-1]) + 1e-6

        # 標準区間 x∈[-1,1] から物理座標 ỹ=L_half(1-x)/2 への写像のヤコビアン。
        # dy/dx = -L_half/2（定数）なので、連鎖律で D_phys = D_std/(dy/dx) と単純にスケールできる。
        dy_dx   = -(L_half / 2.0)
        D_phys  = D_std  / dy_dx
        D2_phys = D2_std / dy_dx**2

        # Boltzmann因子の指数。電子 q_e=-1 → arg_e=+(Φ-β_e A)/Te、
        # イオン q=+1 → arg=-(Φ-β_s A)/T_s と符号が反転する
        # （n_s ∝ exp[-q_s(Φ-β_s A)/T_s] という統一式の帰結）。
        # np.clip で exp のオーバーフローを防ぐ（低温・大振幅で指数の引数が大きくなり得るため）。
        arg_e   = np.clip( (Phi - beta_e   * A) / Te,    -100, 100)
        arg_inc = np.clip(-(Phi - beta_inc  * A) / Tinc,  -100, 100)
        arg_ref = np.clip(-(Phi - beta_ref  * A) / Tref,  -100, 100)

        # 密度の前因子は化学ポテンシャル正規化: 空間平均が背景密度に一致するよう
        # exp[...] をその空間平均で割る。これにより周期境界条件が要求する
        # 積分中性条件（平均電荷中性・平均電流中性）が自動的に満たされる。
        # （比例係数を単純に (1,1-η,η) に固定すると、Neumann境界条件と矛盾し
        #   折り返し点にキンクが生じる。詳細は docs/snr_equilibrium_solver.html）
        e_e = np.exp(arg_e)
        e_i = np.exp(arg_inc)
        e_r = np.exp(arg_ref)
        n_e   = e_e / _avg(e_e)
        n_inc = (1.0 - eta) * e_i / _avg(e_i)
        n_ref = eta          * e_r / _avg(e_r)

        # ガウス則（正規化）: d²Phi/dy² = n_e - n_inc - n_ref
        rhs_Phi = n_e - n_inc - n_ref
        # アンペール則（正規化）: d²A/dy² = beta_e*n_e - beta_inc*n_inc - beta_ref*n_ref
        rhs_A   = beta_e * n_e - beta_inc * n_inc - beta_ref * n_ref

        eq_Phi = D2_phys @ Phi - rhs_Phi
        eq_A   = D2_phys @ A   - rhs_A

        # 全残差に (dy/dx)^2 を掛けて Chebyshev空間相当のスケールに揃える
        # （D2_phys = D2_std/(dy/dx)^2 なので、乗じると (dy/dx)^2 が相殺され
        #  D2_std@A - rhs_A*(dy/dx)^2 の形になる）。L_half 自体が未知数で
        # 反復ごとに変わるため、この再スケーリングで残差各成分の大きさを
        # そろえ、Levenberg-Marquardt の収束を安定させる。
        factor = dy_dx**2
        res = np.zeros_like(u)
        # 内部ノードの方程式（ノード 1..N-1）
        res[1 : N]           = eq_A[1 : N]   * factor
        res[N + 2 : 2*N + 1] = eq_Phi[1 : N] * factor
        # 境界条件（下記コメントの物理的意味は snr_equilibrium_solver.html の図参照）
        res[0]       = (A[0] - current_a0)       * factor   # A(0) = a0
        res[N]       = (D_phys @ A)[N]           * factor   # dA/dy = 0 at y=L_half
        res[N + 1]   = (D_phys @ Phi)[0]         * factor   # dPhi/dy = 0 at y=0
        res[2*N + 1] = (D_phys @ Phi)[N]         * factor   # dPhi/dy = 0 at y=L_half
        res[2*N + 2] = (D_phys @ A)[0]           * factor   # dA/dy = 0 at y=0
        return res

    if u_init is not None:
        # --- warm-start モード ---
        # 既知の近傍解（例: 別のフレーム・別パラメータでの解）を初期値として与え、
        # 線形分散 seed も a0 ランプもスキップして a0_target で1回だけ Newton を解く。
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

        # 小振幅極限 (Φ,A)→0 で密度を1次まで展開すると、次の線形2階ODE系になる:
        #   d²Φ/dy² =  C_PP Φ - C_AP A
        #   d²A/dy² =  C_AP Φ + C_AA A
        # （導出: arg_e≈(Φ-β_e A)/Te 等をexp[...]≈1+argで線形化し、rhs_Phi/rhs_A に代入。
        #   0次（背景一様）成分は電荷・電流中性条件で厳密にキャンセルし、1次成分だけが残る）
        C_PP = 1.0/Te + (1.0-eta)/Tinc + eta/Tref
        C_AP = be_s/Te + (1.0-eta)*bi_s/Tinc + eta*br_s/Tref
        C_AA = -(be_s**2/Te + (1.0-eta)*bi_s**2/Tinc + eta*br_s**2/Tref)

        # 正弦波モード A=cos(ky), Φ=r・cos(ky) を上の線形系に代入すると:
        #   -k²r = C_PP r - C_AP        …(I)  ← Φ方程式から
        #   -k²  = C_AP r + C_AA        …(II) ← A方程式から (k² = -(C_AP r + C_AA))
        # (II) を (I) へ代入して k² を消去すると、比 r=Φ/A が満たす2次方程式が出る:
        #   C_AP r² + (C_AA - C_PP) r + C_AP = 0
        # phi_ratio r が満たす方程式: C_AP*r^2 + (C_AA - C_PP)*r + C_AP = 0
        # （既存コードの a_q=C_APhi, b_q=C_AA-C_PhiPhi, c_q=C_APhi と同じ構造）
        a_q, b_q, c_q = C_AP, C_AA - C_PP, C_AP
        D_q = b_q**2 - 4.0 * a_q * c_q
        if D_q < 0 or abs(C_AP) < 1e-14:
            # 縮退ケース（純静電・ドリフトなし、または C_AP≈0）。
            # r=0（AとΦが結合しない）とし、波数²はC_PP/-C_AAのうち正の方を採用。
            phi_ratio = 0.0
            k2 = -C_AA if -C_AA > 0 else C_PP
        else:
            r1 = (-b_q + np.sqrt(D_q)) / (2.0 * a_q)
            r2 = (-b_q - np.sqrt(D_q)) / (2.0 * a_q)
            # (II) から k²=-(C_AP r+C_AA) を計算し、正（=振動解が存在）の根を採用。
            phi_ratio, k2 = r1, -(C_AP * r1 + C_AA)
            if k2 <= 0:
                phi_ratio, k2 = r2, -(C_AP * r2 + C_AA)
            if k2 <= 0:
                # 両根とも k²≤0 なら、指数的に発散/減衰するだけで振動する
                # （=周期的フィラメント解が存在する）根がない。
                raise ValueError(
                    "指定パラメータで振動解が存在しません。"
                    f"C_PP={C_PP:.4f}, C_AP={C_AP:.4f}, C_AA={C_AA:.4f}"
                )

        # 半周期の初期推定値: A=cos(ky) が y=0 で最大、y=L_half で dA/dy=0 となる
        # （境界のNeumann条件を満たす）ためには k・L_half=π、すなわち L_half=π/k。
        L_half_guess = np.pi / np.sqrt(k2)
        # cos_profile(x) = a0_start・cos(π(1-x)/2) は、写像 y=L_half(1-x)/2 の下で
        # a0_start・cos(π y/L_half) = a0_start・cos(k y) に等しい（k=π/L_half を代入）。
        # つまり x_std 上でこの式を評価するだけで、物理空間の cos(ky) モードが作れる。
        cos_profile = a0_start * np.cos(np.pi * (1.0 - x_std) / 2.0)
        u_guess = np.zeros(2 * N + 3)
        u_guess[:N + 1]          = cos_profile        # A の初期推定（振幅 a0_start）
        u_guess[N + 1 : 2*N + 2] = phi_ratio * cos_profile  # Phi = r・A の初期推定
        u_guess[-1]              = L_half_guess

        # --- Newton継続法: a0 を a0_start から a0_target まで段階的に増大 ---
        # 大振幅では直接 Newton を投げても収束しないことが多いため、
        # 小振幅の線形解から出発し、a0 を刻みながら「前ステップの解」を
        # 次ステップの初期値にする継続法（homotopy）で頑健に収束させる。
        a0_ramp = np.linspace(a0_start, a0_target, n_steps)
    for i, current_a0 in enumerate(a0_ramp):
        if i > 0:
            # 線形応答では振幅がa0に比例するはずなので、前回の解を
            # 「今回のa0/前回のa0」倍にスケールしてから Newton に投げると
            # 初期残差が小さく収束が速い（L_half は振幅に依らずそのまま流用）。
            u_guess[:2*N + 2] *= current_a0 / a0_ramp[i - 1]
        # scipy.optimize.root の Levenberg-Marquardt法。過決定・特異点近傍でも
        # 比較的頑健に動くため、非線形スペクトル残差の求解に採用。
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

    # 場の定義（符号規約は snr_linear.py と共通）。
    B0z = -(D_phys @ A_sol)    # B_{0z} = -dA/dy
    E0y = -(D_phys @ Phi_sol)  # E_{0y} = -dPhi/dy

    # 物理y座標（半周期、c/omega_pe単位）。x_std は1(y=0)から-1(y=L_half)へ並ぶので
    # y_half = L_half(1-x)/2 は 0 から L_half へ単調増加する。
    y_half = L_half * (1.0 - x_std) / 2.0

    # 密度プロファイル（残差と同じ化学ポテンシャル正規化を適用）
    arg_e   =  (Phi_sol - beta_e   * A_sol) / Te
    arg_inc = -(Phi_sol - beta_inc  * A_sol) / Tinc
    arg_ref = -(Phi_sol - beta_ref  * A_sol) / Tref
    e_e = np.exp(arg_e)
    e_i = np.exp(arg_inc)
    e_r = np.exp(arg_ref)
    n_e_half   = e_e / _avg(e_e)
    n_inc_half = (1.0 - eta) * e_i / _avg(e_i)
    n_ref_half = eta          * e_r / _avg(e_r)

    # 1周期分に対称拡張:
    #   B, E: 交互フィラメントのため反対称（奇関数）
    #   密度: 対称（偶関数）
    # v[-2::-1] は「最後の要素を除いて逆順」＝ 折り返し点(索引 N)を二重にしないための
    # スライス（v=[v0,...,vN] なら v[-2::-1]=[v(N-1),...,v0]）。
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
        "lambda0":   2.0 * L_half,
        "beta_e":    beta_e,
        "u_solution": sol.x.copy(),
        "params": {
            "eta": eta, "beta_inc": beta_inc, "beta_ref": beta_ref,
            "Te": Te, "Tinc": Tinc, "Tref": Tref, "a0_target": a0_target,
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
) -> dict:
    """電子静止系（β_e=0）の周期平衡を最初から直接解く。

    支配方程式は solve_snr_equilibrium と同一だが、電子ドリフトを 0 に固定して
    解く。boost や frame-homotopy（自然系で解いてから系を移す2段階）は行わず、
    a0 ランプ1本で電子系平衡を直接求める。

    入力の (beta_inc, beta_ref) は自然系（電流中性系）のドリフトを渡す。内部で
    β_e,nat = (1-eta)*beta_inc + eta*beta_ref を差し引いた電子系ドリフト
    (beta_inc-β_e,nat, beta_ref-β_e,nat) に移す。これらは (1-eta)β_inc'+η β_ref'=0
    を満たし、solve_snr_equilibrium が再導出する beta_e も 0 になる
    （ガリレイ的な速度シフトなので電流中性条件はシフト後も自動的に保たれる:
      (1-η)(β_inc-β_e,nat)+η(β_ref-β_e,nat) = [(1-η)β_inc+ηβ_ref] - β_e,nat = 0）。

    ただし β_e=0 では線形分散 seed の結合係数 C_AP が縮退（Tinc=Tref で C_AP=0、
    なぜなら電子系ドリフトはすでに電流中性 (1-η)bi_e+η br_e=0 を満たすため、
    Tinc=Tref のとき C_AP=(bi_e(1-η)+br_e η)/T=0 になる）し λ0 が崩壊した
    spurious 解に落ちる。これを避けるため seed のみ自然系ドリフト
    (_seed_drifts) で初期化して正しい basin を選ぶ。残差・抽出（E0y など）は
    β_e=0・電子系ドリフトで行うので、返る E0y は電子系の物理場そのもの。

    Parameters
    ----------
    eta, beta_inc, beta_ref, Te, Tinc, Tref, a0_target, N_points, n_steps:
        solve_snr_equilibrium と同じ。beta_inc/beta_ref は自然系の値。

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
        _seed_drifts=(beta_inc, beta_ref, beta_e_nat),
    )
