# マッハ数基準のパラメータ構築

代表パラメータとしてマッハ数を1つ決めれば、速度 → 温度 → 各種 β の順で
3成分SNRプラズマ平衡のパラメータが一意に決まる。本ノートはその導出チェーンと、
現行 PIC ラン (`config_2d_inp_xy.toml`) での検算をまとめる。実装は
`fmi.mach_parameters`（`MachConfig`, `build_solver_params`, `build_pic_config`）。

## 1. 記号と正規化

| 記号 | 意味 |
|------|------|
| $M_A$ | イオンAlfvénマッハ数 $\beta_{sh}/\beta_{A,i}$ |
| $M_S$ | イオン音速マッハ数 $\beta_{sh}/\beta_{cs,i}$ |
| $\eta$ | 反射イオン割合（PIC の `alpha`） |
| $m_i/m_e$ | 質量比（`mime`） |
| $\sigma$ | 磁化 $(v_{A,e}/c)^2 = (\Omega_{ce}/\omega_{pe})^2$。$\sigma=0$ で非磁化 |
| $\gamma$ | 断熱指数（既定 $5/3$） |

**速度記法（統一規約）:** すべての速度は $c$ で正規化し $\beta_{\langle名\rangle}\equiv v/c$ と綴る。
ドリフト $\beta_e,\beta_{inc},\beta_{ref}$、熱速度 $\beta_{th,s}$、音速 $\beta_{cs,s}=\sqrt\gamma\,\beta_{th,s}$、
衝撃波 $\beta_{sh}$、Alfvén $\beta_{A,s}$ すべて同じ規約。マッハ数は速度比なので $c$ で割っても同値
（$V_{sh}/v_{A,i}=\beta_{sh}/\beta_{A,i}$）。
**別物:** 「プラズマ $\beta$」（圧力比 $2nk_BT/(B^2/2\mu_0)$）は速度ではない。コードでは圧力比のみ
`beta_e_plasma` と綴る。

正規化は `snr_equilibrium` と同一：長さ $c/\omega_{pe}$、温度 $m_e c^2/k_B$、速度 $\beta=v/c$。
電子質量 $m_e=1$、イオン質量 $m_i=$ `mime`。
規格化温度は $\tilde T_s = k_B T_s/(m_e c^2)$ で、熱速度と $\tilde T_{inc}=m_i\,\beta_{th,i}^2$、
$\tilde T_e=\beta_{th,e}^2$ で結ばれる。

## 2. 導出チェーン（磁化時）

$$
\begin{aligned}
\text{(1)}\quad & \beta_{A,i} = \sqrt{\sigma/(m_i/m_e)} \\
\text{(2)}\quad & \beta_{sh} = M_A \cdot \beta_{A,i} \\
\text{(3)}\quad & \beta_{d} = \beta_{sh}/2 \qquad(\text{相対ドリフトの半分 } 2\beta_{d}=\beta_{sh}) \\
\text{(4)}\quad & \beta_{inc} = \frac{-2\eta\,\beta_{d}}{1-(1-2\eta)\beta_{d}^2},\quad
                   \beta_{ref} = \frac{+2(1-\eta)\,\beta_{d}}{1+(1-2\eta)\beta_{d}^2} \\
\text{(5)}\quad & \beta_{cs,i} = \beta_{sh}/M_S \\
\text{(6)}\quad & \tilde T_{inc} = (m_i/m_e)\,\beta_{cs,i}^2/\gamma
\end{aligned}
$$

(4) は鏡面反射規約のシム系ドリフトで、PIC `main.cpp:61-62`（`vdi, vdr`）と同式。
電子ドリフトはソルバー内で電流中性条件
$\beta_e=(1-\eta)\beta_{inc}+\eta\beta_{ref}$ から自動導出される（$\approx 0$）。

電子温度は **独立ノブ**：`beta_e_plasma`（プラズマβ）/ `beta_te`（$\beta_{th,e}$）/ $T_e/T_{inc}$ の
いずれか1つを指定する。プラズマβ 指定時は
$\tilde T_e = \beta_{th,e}^2 = \sigma\,\beta_e^{plasma}/2$。

## 3. マッハ恒等式

プラズマ $\beta$（圧力比）の定義 $\beta_s^{plasma} = 2 \beta_{th,s}^2/\beta_{A,s}^2$ と
音速 $\beta_{cs}=\sqrt{\gamma}\,\beta_{th}$ から、

$$
\frac{M_A}{M_S} = \frac{\beta_{cs,i}}{\beta_{A,i}} = \sqrt{\frac{\gamma\beta_i}{2}}
\quad\Longrightarrow\quad
\beta_i = \frac{2}{\gamma}\left(\frac{M_A}{M_S}\right)^2 .
$$

これが指導教員の式 $\beta=(2/\gamma)(M_A/M_S)^2\sim(M_A/M_S)^2$ に対応する。

## 4. 磁化 / 非磁化の分岐

- **磁化 ($\sigma>0$)**: $M_A$ が速度を、$M_S$ が温度をアンカーする。
  両方が意味を持ち、$\beta_i$ は上の恒等式で決まる。
- **非磁化 ($\sigma=0$)**: $\beta_A\to 0$ で $M_A\to\infty$ となり無意味。
  $\beta_{sh}=V_{sh}/c$ を直接アンカー（コード `beta_sh`）し、$M_S$ のみが温度を決める
  （「背景磁場を考えないなら音速マッハ数のみ考えれば良い」に対応）。
  `build_pic_config` は `needs_direct_temperature=True` を返し、温度を $\beta_A$ から
  切り離して直接指定する PIC 初期化（`main.cpp` 小改修）が必要であることを示す。

## 5. 検算表（現行ラン）

`config_2d_inp_xy.toml`: $\sigma=0.0025$, $m_i/m_e=400$, `ush`$=0.125$,
$\eta=\alpha=0.2$, `betae`$=8$, `betai`$=$`betar`$=0.5$, $\gamma=5/3$。
これは $M_A=100$, $M_S\approx 155$（厳密 $154.92$）に対応する。

| 量 | `build_solver_params` 出力 | config 実値 |
|----|------|------|
| $M_A$ | 100（入力） | 100（コメント記載） |
| $\beta_{A,i}$ | 0.0025 | $\sqrt{0.0025/400}=0.0025$ |
| $\beta_{sh}$ | 0.25 | $2\beta_{d}\approx0.25$ |
| $\beta_{inc}$ | $-0.0505$ | — |
| $\beta_{ref}$ | $+0.1981$ | — |
| $\beta_e$ | $-7.5\times10^{-4}\approx0$ | 電子静止系で $\approx0$ |
| $\tilde T_{inc}$ | $6.243\times10^{-4}$ | `betai`$=0.5\Rightarrow6.25\times10^{-4}$ |
| $\tilde T_e$ | $0.01$ | `betae`$=8\Rightarrow \beta_{th,e}=0.1$ |

`build_pic_config` の逆生成: `betae`$=8.0$, `betai`$\approx0.50$, `ush`$\approx0.126$
（$M_S=155$ の丸めで `betai`/$\tilde T_{inc}$ は約 0.1% ずれる。厳密 $M_S=154.92$ で一致）。

回帰テスト: `tests/test_mach_parameters.py`。
