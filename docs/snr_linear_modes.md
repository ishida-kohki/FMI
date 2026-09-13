# SNR平衡の線形固有値解析（非相対論 warm 多流体）

SNR 初期平衡（電子・入射イオン・反射イオン）を背景に、warm 多流体 + Maxwell を
線形化し、多周期にタイルした領域で Fourier スペクトル直接固有値解析を行う。目的は
FMI（フィラメント併合）の成長モードを示すこと。正規化は `snr_equilibrium` と同一
（長さ c/ω_pe、速度 β=v/c、c=1、m_e=1, m_i=mime、ω_pe=1）。

## 1. 幾何と偏波分離

平衡は y のみに依存。電流は x（ドリフト方向）、磁場 B0z は面外 (z)。摂動は

    δg(x, y, t) = δĝ(y) · exp(i kx x − i ω t)

kx はパラメータ、y はスペクトルで解像、z 方向一様。x-y 面では2つの偏波が分離する。
不安定（FMI/DKI）に効くのは **P 偏波** {δE_x, δE_y, δB_z} と流体 {δn_s, δv_sx, δv_sy}。
O 偏波 {δE_z, δB_x, δB_y, δv_sz} は分離するので除外する。

## 2. 背景量

各成分 s ∈ {e, inc, ref}:
- 密度 n_{s0}(y)（y で変動、平衡ソルバーの ne/ninc/nref）
- ドリフト u_{s0x} = β_s（**y 方向に一定**）
- 温度 T_s（一定）、圧力 p_{s0} = n_{s0} T_s、摂動圧力 δp_s = γ_ad T_s δn_s
- 質量 m_e=1, m_inc=m_ref=mime、電荷 q_e=−1, q_inc=q_ref=+1

**u_{s0x} が y 一定であることの確認:** Boltzmann 閉包 n_{s0} ∝ exp(∓(Φ − β_s A)/T_s)
より T_s d(ln n_{s0})/dy = E0y − β_s B0z。一方 y 方向の力釣り合いは
q_s(E0y − u_{s0x} B0z) = (1/n_{s0}) dp_{s0}/dy = T_s d(ln n_{s0})/dy。
両者が一致するのは u_{s0x} = β_s c（一定）のとき。

**背景場は線形摂動式にも現れる（2026-07 訂正）:** 密度揺らぎ δn_s の分だけ粒子が
増減すると、その分も背景場から力を受ける。この一次の力 q_s δn_s (E0y − β_s B0z)
（力釣り合いより T_s δn_s d(ln n_{s0})/dy に等しい）が運動 y に入る。
旧版はこの項を欠いており、「E0y は線形摂動式に直接現れない」としていたのは誤り。
欠落は並進ゼロモード検査（下記 §5-3）で検出・確認した。

**フレーム（電流中性系を主に使う）:** β_e は独立量ではなく電流中性条件
β_e = (1−η)β_inc + η β_ref でドリフトから決まる（既定で β_e ≈ −7.5e-4）。フレーム変換は
ガリレイ・ブースト Φ→Φ−β_e A, β_s→β_s−β_e で、密度・B0z・相対電流は不変、E0y のみ
E0y→E0y−β_e B0z だけずれる。背景場力項の係数 (E0y − β_s B0z) はこの変換で不変
（(E0y−β_eB0z) − (β_s−β_e)B0z = E0y − β_sB0z）。ドリフトの一様シフト β_e は
Doppler 因子 Ω_s = ω − kx β_s を ω の一定シフト kx β_e に移すだけなので、線形スペクトルは
フレーム共変（β_e~1e-4 で無視可）。よって背景は電流中性系（solve_snr_equilibrium）で供給すれば十分で、
電子静止系（β_e=0）は PIC フレーム合わせ用の付録 boost にすぎない。

## 3. 線形化方程式（12場、∂_t → −iω で ω 一次）

Doppler 因子 Ω_s ≡ ω − kx u_{s0x}。c=1。

連続:
    −iω δn_s + i kx(n_{s0} δv_sx + u_{s0x} δn_s) + ∂_y(n_{s0} δv_sy) = 0

運動 x:
    m_s n_{s0}(−iΩ_s) δv_sx = q_s n_{s0} δE_x + q_s n_{s0} B0z δv_sy − i kx γ_ad T_s δn_s

運動 y:
    m_s n_{s0}(−iΩ_s) δv_sy = q_s n_{s0} δE_y − q_s n_{s0}(B0z δv_sx + u_{s0x} δB_z)
                              + q_s δn_s (E0y − β_s B0z) − ∂_y(γ_ad T_s δn_s)

Faraday:
    −iω δB_z = −(i kx δE_y − ∂_y δE_x)

Ampère x:
    −iω δE_x = ∂_y δB_z − Σ_s q_s(n_{s0} δv_sx + u_{s0x} δn_s)

Ampère y:
    −iω δE_y = −i kx δB_z − Σ_s q_s n_{s0} δv_sy

電流結合係数は ω_pe=1 正規化で 1（Gaussian の 4π を吸収）。∇·E=ρ（Poisson）は
連続 + Ampère から保たれるので検証に使う。

## 4. 標準固有値問題への整形と離散化

各式を ∂_t 項の係数で割ると ω x = L x（標準固有値）になる。係数:
- 連続: −i、運動: −i m_s n_{s0}、Faraday/Ampère: −i。
1/(−i)=i を使い、L の各ブロックは下記（ω·場 = Σ ブロック · 場）:

連続 (δn_s):   δn_s: kx u_{s0x},  δv_sx: kx n_{s0},  δv_sy: −i ∂_y(n_{s0}·)
運動x (δv_sx): δv_sx: kx u_{s0x},  δE_x: i q_s/m_s,  δv_sy: i q_s B0z/m_s,
               δn_s: (kx γ_ad T_s/m_s)(1/n_{s0})
運動y (δv_sy): δv_sy: kx u_{s0x},  δE_y: i q_s/m_s,  δv_sx: −i q_s B0z/m_s,
               δB_z: −i q_s u_{s0x}/m_s,
               δn_s: −i(γ_ad T_s/m_s)(1/n_{s0})∂_y + i(q_s/m_s)(E0y − β_s B0z)/n_{s0}
Faraday (δB_z): δE_y: kx,  δE_x: i ∂_y
Ampère x (δE_x): δB_z: i ∂_y,  δv_sx: −i q_s n_{s0},  δn_s: −i q_s u_{s0x}
Ampère y (δE_y): δB_z: kx,  δv_sy: −i q_s n_{s0}

**離散化:** [0, Nλ0] 上の一様 M 点で Fourier スペクトル。∂_y → Fourier 微分行列 D_y。
背景 n_{s0}(y), B0z(y), E0y(y) は実空間で対角（np.diag）。E0y は電子の力釣り合い
E0y = β_e B0z − T_e n_e'/n_e から格子上で再構成する（Boltzmann 平衡で厳密）。
各場を M ベクトルとして (3·成分数 + 3) ブロック、総サイズ (3S+3)·M。
`numpy.linalg.eig`（密、全スペクトル）で解き、Im(ω)>0 を成長モードとして抽出、
解像度収束で物理固有値を選別。

**密度下限（低温で必須）:** 運動 y の 1/n_{s0}（圧力項・背景場力項）は密度が
空乏化する低温平衡（n_ref → ~1e-54）で発散する。`density_floor` により
1/max(n_{s0}, floor·max n_{s0}) に置換する。斜め枝は floor 2桁変化に対し γ 変動
~5%（系統誤差として報告）。kx=0 の FMI は floor 依存が強く、この流体閉包では
頑健に決まらない（空乏域では Boltzmann 流体近似自体が適用限界）。

## 5. 検証方針

1. **光波分散較正（一様・単一電子）:** ω² = ω_pe² + c²kx²（= 1 + kx²）を再現 → EM+
   流体結合と結合係数が正しい。
2. **cold filamentation（一様・逆流2ビーム）:** 解析分散
   γ²(k_y) = [√((k_y²+ω_pe²)² + 4β²k_y²ω_pe²) − (k_y²+ω_pe²)] / 2
   （ω_pe=1）を再現 → ドリフト・磁気結合（v0×δB → δn → 電流）のループが正しい。
   k_y→∞ で γ → β ω_pe。
3. **並進ゼロモード検査（非一様背景の較正）:** 平衡の y 微小並進
   （δn_s=−n0_s', δEy=−E0y', δBz=−B0z', δv=0）は kx=0, ω=0 の厳密解であり、
   演算子はこれを消すべき。1・2 の較正は一様背景のみで勾配結合項を検査できない
   ため、この検査が背景場力項・圧力項の勾配部分を較正する
   （tests/test_snr_linear.py::test_translation_zero_mode_annihilated）。
4. **SNR-FMI:** SNR 平衡背景・kx=0 で Im(ω)>0 の成長モードが存在し、最速モードの
   δB_z が subharmonic（k_y/k0 < 1 = 隣接フィラメント併合）。
5. **収束:** 成長率が points_per_period, n_periods に対し収束（spurious 除外）。
6. **Poisson 残差:** 固有ベクトルで ∇·E − ρ が小さいこと（任意）。
