# SGEPSS 2026 予稿

- 学会: SGEPSS 第160回総会・講演会（2026年11月9-12日）
- 予稿締切: 2026年8月19日
- 提出形式: 日本語・英語の両方（各2000字）
- 最終更新: 2026年8月14日

---

## 日本語版

　無衝突衝撃波は長年にわたり研究されてきた。その主な動機は、荷電粒子の加速機構の理解である。中でも、衝撃波遷移層におけるイオンビーム不安定性は、電子加速問題を解く鍵として重要視されてきた。超新星残骸 (SNR) のような高マッハ数領域では、Weibel instability (WI) が支配的なモードとして注目されている [Medvedev & Loeb, 1999; Kato & Takabe, 2008; Matsumoto et al., 2017]。衝撃波遷移層では、WIに伴い、電流層が自発的に生成し、磁場の増幅が起こる。その後、電流層の崩壊にともなって電子の加熱・加速が起こると考えられている [Matsumoto et al., 2015]。

　WIが衝撃波遷移層に電流層を形成することはすでに知られているが、その後どのように電流層が崩壊するかは明らかではない。SNR衝撃波を模したPICシミュレーションでは、これらの電流層上で自発的な磁気リコネクションが生じることが示されている [Matsumoto et al., 2015; Jikei et al., 2024]。また、一般に電流層の不安定性として磁気リコネクションを起こす Tearing Instability (TI) に加えて Drift Kink Instability (DKI) [Pritchett et al., 1996] が重要な役割を果たすことが知られており、Drift Sausage Instability (DSI) の存在も指摘されている [Yoon & Lui, 2001]。また、WIの線形段階で形成する電流フィラメントは Filament Merging Instability (FMI) [Vanthieghem et al., 2018] で長波長に進化することも知られている。このように電流層崩壊を引き起こしうる不安定性は複数あるが、その役割はよく分かっていない。

　電流層に対する既存の線形解析の多くは、実際にSNR衝撃波遷移領域に形成される自己無撞着な状態ではなく、Harris sheet のような理想化された平衡に依拠している。そのため、SNR衝撃波において、これらの不安定性のうちどれが実際に生じるのかは明らかでない。Vanthieghem et al. [2018] は、Weibel不安定性の結果として生じる電流フィラメントの非線形発展の理解を目的として、相対論的ペアプラズマにおけるDKIとFMIの競合を調べた。本研究では、これをSNR衝撃波に対応する非相対論的なイオン・電子系で行うことを目指す。

　本研究は2段階で進める。第一に、電子・入射イオン・反射イオンの3流体モデルを用いて、周期的な平衡解を数値的に構築する。これはWIがある程度成長した結果として生じる準平衡状態をモデル化したものである。第二に、この平衡のまわりで線形固有値問題を解き、TI、DKI、FMI、DSI などの各種不安定性の成長率のパラメータ依存性を調べ、WIの非線形発展段階における電流層の安定性を議論する。

　これまでに特定のパラメータ設定で平衡解を数値的に構築し、線形解析を行った。その結果、波数ベクトルが電流に平行な場合、DKI と DSI の二つの不安定性が現れることがわかった。

　今後の課題は複数ある。まず、線形解析のパラメータ依存性を明らかにしたい。また、TIやFMIが支配的になる可能性を探るため、電流に平行でない波数ベクトルへの拡張も必要である。さらに、シミュレーションを行って、線形解析の結果を検証し、不安定性の非線形発展を調べる予定である。

---

## English version

　　Collisionless shocks have been studied for decades. The primary motivation is to understand how charged particles are accelerated. Among the relevant processes, ion beam instabilities in the shock transition region have been regarded as the key to the electron acceleration problem. In the high Mach number regime relevant to supernova remnants (SNRs), the Weibel instability (WI) has attracted attention as the dominant mode [Medvedev & Loeb, 1999; Kato & Takabe, 2008; Matsumoto et al., 2017]. In the shock transition region, WI spontaneously forms current sheets and amplifies the magnetic field. Electrons are then thought to be heated and accelerated as these current sheets break up [Matsumoto et al., 2015].

　　It is already known that WI forms current sheets in the shock transition region, but how these sheets subsequently break up remains unclear. PIC simulations modeling SNR shocks have shown that spontaneous magnetic reconnection occurs at these current sheets [Matsumoto et al., 2015; Jikei et al., 2024]. More generally, the tearing instability (TI), which drives magnetic reconnection, and the drift kink instability (DKI) [Pritchett et al., 1996] are known to play important roles as current sheet instabilities, and the existence of a drift sausage instability (DSI) has also been pointed out [Yoon & Lui, 2001]. In addition, the current filaments formed during the linear stage of WI are known to evolve toward longer wavelengths through the filament merging instability (FMI) [Vanthieghem et al., 2018]. Thus, several instabilities can drive the breakup of the current sheets, but their respective roles are not well understood.

　　Most existing linear analyses of current sheets rely on idealized equilibria such as the Harris sheet, rather than on the self-consistent state that actually forms in the SNR shock transition region. It therefore remains unclear which of these instabilities operates in SNR shocks. Vanthieghem et al. [2018] examined the competition between DKI and FMI in relativistic pair plasmas, aiming to understand the nonlinear evolution of the current filaments produced by WI. This study aims to carry out the corresponding analysis for the non-relativistic ion-electron system relevant to SNR shocks.

　　We proceed in two stages. First, we numerically construct a periodic equilibrium using a three-fluid model of electrons, incoming ions, and reflected ions. This models the quasi-equilibrium state that results once WI has grown to a certain amplitude. Second, we solve the linear eigenvalue problem about this equilibrium, examine how the growth rates of the individual instabilities, namely TI, DKI, FMI, and DSI, depend on the parameters, and discuss the stability of the current sheets in the nonlinear stage of WI.

　　So far we have numerically constructed the equilibrium for a particular parameter set and performed the linear analysis. We find that two instabilities appear when the wavevector is parallel to the current: DKI and DSI.

　　Several tasks remain. First, we want to clarify the parameter dependence of the linear analysis. Second, an extension to wavevectors not parallel to the current is necessary in order to explore the possibility that TI or FMI becomes dominant. Finally, we plan to run simulations to verify the results of the linear analysis and to investigate the nonlinear evolution of the instabilities.

---

## 引用文献

すべて実物PDF、arXiv、または出版社ページで書誌を確認済み。

```
Medvedev, M. V., & Loeb, A. (1999). Generation of magnetic fields in the
relativistic shock of gamma-ray burst sources.
The Astrophysical Journal, 526, 697-706.

Kato, T. N., & Takabe, H. (2008). Non-relativistic collisionless shocks in
unmagnetized electron-ion plasmas.
The Astrophysical Journal Letters, 681, L93-L96. arXiv:0804.0052

Matsumoto, Y., Amano, T., Kato, T. N., & Hoshino, M. (2017). Electron surfing
and drift accelerations in a Weibel-dominated high-Mach-number shock.
Physical Review Letters, 119, 105101. arXiv:1709.03673

Matsumoto, Y., Amano, T., Kato, T. N., & Hoshino, M. (2015). Stochastic electron
acceleration during spontaneous turbulent reconnection in a strong shock wave.
Science, 347(6225), 974-978.

Jikei, T., Amano, T., & Matsumoto, Y. (2024). Enhanced magnetic field
amplification by ion-beam Weibel instability in weakly magnetized
astrophysical shocks.
The Astrophysical Journal, 961, 157. doi:10.3847/1538-4357/ad1594

Pritchett, P. L., Coroniti, F. V., & Decyk, V. K. (1996). Three-dimensional
stability of thin quasi-neutral current sheets.
Journal of Geophysical Research, 101(A12), 27413-27429. doi:10.1029/96JA02665

Yoon, P. H., & Lui, A. T. Y. (2001). On the drift-sausage mode in
one-dimensional current sheet.
Journal of Geophysical Research, 106(A2), 1939-1947. doi:10.1029/2000JA000130

Vanthieghem, A., Lemoine, M., & Gremillet, L. (2018). Stability analysis of a
periodic system of relativistic current filaments.
Physics of Plasmas, 25, 072115. doi:10.1063/1.5033562
```

### 引用の根拠

| 引用 | 支える記述 | 確認した文言 |
|---|---|---|
| Medvedev & Loeb 1999 | 高マッハ数WIの起点 | GRB相対論的衝撃波が対象。Jikei 2024 も相対論的研究群の先頭として引用 |
| Kato & Takabe 2008 | 非相対論・電子イオン系への拡張、SNRへの適用 | "Weibel-mediated collisionless shocks are driven at non-relativistic propagation speed (0.1c < V < 0.45c) in unmagnetized electron-ion plasmas... probably explains the robust formation of collisionless shocks, for example, driven by young supernova remnants" |
| Matsumoto et al. 2017 | 高マッハ数でWIが支配的 | 題名が "Weibel-dominated High-Mach-number Shock" |
| Matsumoto et al. 2015 | 電流層上の自発的磁気リコネクション、電子加熱・加速 | "supercomputer simulations showing..." (Science 347, 974) |
| Jikei et al. 2024 | SNR衝撃波を模したPICでの電流層形成とリコネクション | "We conclude that this scenario is applicable to typical young supernova remnant shocks" |
| Pritchett et al. 1996 | DKI（およびTI） | 題名 "Three-dimensional stability of thin quasi-neutral current sheets"。DKIのパリティ構造とtearing的モードの両方を扱う |
| Yoon & Lui 2001 | 電流層のsausageモードの存在 | "the charge neutrality condition leads to sausage-type fluctuations, which propagate along the direction of the cross-field current flow, hence the drift-sausage mode" |
| Vanthieghem et al. 2018 | FMI、および相対論的ペアプラズマでのDKI-FMI競合 | 本研究が直接対応づける先行研究 |

---

## 字数

### 日本語版

制限内。

### 英語版

| 段落 | 字数（空白込） | 語数 |
|---|---|---|
| 1 | 714 | 109 |
| 2 | 920 | 141 |
| 3 | 593 | 89 |
| 4 | 537 | 83 |
| 5 | 239 | 37 |
| 6 | 388 | 62 |
| **合計** | **3391** | **521** |

**2000字制限を1391字超過。未対応。**

削減候補:

1. 第2段落1文目が第1段落末と重複（"It is already known that WI forms current sheets in the shock transition region" → 削除して "However, how these sheets subsequently break up remains unclear." に）。約90字減。
2. 第1段落1-2文目を統合。約40字減。
3. 第4段落の不安定性の固有名を落とし "the individual instabilities" のみに。約50字減。

上記を全部やっても約1200字残る。先生が明示的に要望した内容（引用群、2段階の説明、今後の課題の分量）が超過の主因なので、どの段落を落とすかは要相談。

---

## 未解決の指摘

執筆時に検証して指摘したが、本文には反映していない項目。

### 1. TIは既にK=0の解析に含まれているはず

第6段落は「TIやFMIが支配的になる可能性を探るため、電流に平行でない波数ベクトルへの拡張も必要」と書いているが、tearing modeの波数は電流方向（＝反平行磁場方向）に立つ。これは本研究の K = 0 の設定そのものであり、TIは既に今回の解析に含まれているはず。斜め方向で新たに現れるのはFMI（フィラメント合体方向＝電流に垂直）と中間の斜めモード。

### 2. 「DSI」の同定が未確定

上記1の帰結として、K = 0 で得た偶パリティ枝がtearingである可能性を現時点で排除できていない。tearingの磁場摂動は電流層中心について対称。

加えて Vanthieghem et al. 2018 との用語衝突がある。同論文は "sausage-type" を独立の不安定性ではなく副次DKI枝の見え方として扱う。

> "The perturbed magnetic (resp. density) profile of both eigenmodes is even (resp. odd) with respect to the filament center, which is indicative of DKI. ... For the sub-dominant mode, the out-of-phase oscillations of adjacent current filaments translate into sausage-type magnetic fluctuations."

さらに、名称を借りた Yoon & Lui 2001 の drift-sausage mode は **charge quasi-neutrality を仮定した場合にのみ現れる**枝であり、同じ仮定の下では kink 解が存在しない。同論文は結論で "the assumption of quasi-neutrality may be a poor choice in describing the stability of Harris current sheet equilibrium" と自ら述べている。本研究は Poisson を解いており quasi-neutral 仮定を置かないため DKI と偶パリティ枝の両方が出る。名称の借用は妥当だが、Yoon & Lui の drift-sausage mode と同一のモードかは未確認。

安全側に倒すなら第5段落を次のように書き換える:

> その結果、波数ベクトルが電流に平行な場合、フィラメント中心に対して奇パリティ・偶パリティの二つの不安定分枝が現れることがわかった。前者は DKI に対応する。後者は DSI あるいは TI に対応しうるが、その同定は今後の課題である。

### 3. 「電子加速問題」の用語

分野の標準語は「電子入射問題 (electron injection problem)」で、Jikei et al. 2024 も injection を使う。変えるなら日英同時に。

### 4. 「電流層」と「電流フィラメント」の呼び分け

2次元断面ではフィラメント、1次元断面では層で同じものを指す。第2段落に一言（「WIが作る電流フィラメントは、断面で見れば電流層をなす」など）入れると読者に伝わる。

### 5. 「入射イオン」の英訳

`incoming ions` を採用。手元の論文ライブラリ全体で `incident ion` は0件、`incoming ion` は Jikei 2024 に3件・Amano 2022 に1件。Jikei et al. 2024 は同じ3種構成を "the background electrons, the reflected ions, and the incoming ions" と呼ぶ。`incident` は incident shock / incident wave と衝突するため回避。
