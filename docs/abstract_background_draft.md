# SJEPS abstract — 背景パート ドラフト（WI→SNR）

締切 2026-08-19（2000字・英語）。本ファイルは abstract の**背景／動機パート**のみ（方法・結果は別途）。

> **引用の注意（citation-verification）**: 下記引用は well-known か ~/paper/ で確認したもの。abstract 提出前に各文献の正確な書誌・主張を最終確認すること。特に SNR 反射イオン・フィラメント化の代表引用は、ご自身の方針で選定・確定してください。

---

## 案A（5文・標準／推奨, writing-anti-ai 適用済）

> Magnetic-field amplification at collisionless shocks is essential for cosmic-ray acceleration in young supernova remnants (SNRs), but how the required fields are first generated remains unclear. In the collisionless regime, the Weibel (current-filamentation) instability converts the free energy of counter-streaming plasmas into magnetic field, and it has been invoked for both relativistic (e.g., gamma-ray bursts) and non-relativistic shocks [Weibel 1959; Medvedev & Loeb 1999]. At an SNR shock, ions reflected at the front stream back into the incoming plasma; this counter-streaming drives beam filamentation, breaking the precursor into quasi-periodic current filaments. Their growth and merging, and hence the magnetic structure that survives, are set by the linear stability of the filaments, mainly the filamentation/merging (FMI) and drift-kink (DKI) branches. Yet the linear spectrum of a self-consistent SNR filament equilibrium, rather than the idealized Harris sheet used in classical drift-kink studies, has not been characterized.

（→ ここから方法文へ: "Here we construct the three-species SNR current-filament equilibrium and solve its P-polarization linear eigenvalue problem ..."）

## 案B（3文・圧縮／字数が厳しい場合, writing-anti-ai 適用済）

> Magnetic-field amplification at collisionless shocks is central to cosmic-ray acceleration in young supernova remnants (SNRs). The Weibel (current-filamentation) instability converts counter-streaming free energy into magnetic field at both relativistic and non-relativistic shocks [Weibel 1959; Medvedev & Loeb 1999]; at an SNR shock, the counter-streaming of reflected and incoming ions breaks the precursor into current filaments. Their growth and merging, which set the surviving magnetic structure, are governed by the linear stability of the filaments (FMI and drift-kink branches), whose spectrum for a self-consistent, non-Harris SNR equilibrium has not been characterized.

## 推敲メモ（何を直したか）
- filler の adverb "efficiently" を削除。
- promotional な "is now established as a primary field-generation channel" → "has been invoked for"（誇張を避け honest に）。
- rule-of-three "grow, merge, and set the observed magnetic structure" → "growth and merging, and hence the magnetic structure"（三段列挙を解消）。
- em-dash の多用（"—principally"、"—non-Harris—"、"—as opposed to—"）を全廃し comma/appositive に。
- "underpins / seeded and sustained / unsettled" → "essential / first generated / unclear"（平易・直接）。
- リズム: 3文目に semicolon を入れ長さを変化。

---

## 日本語グロス（案Aの意訳）
1. SNR の無衝突衝撃での磁場増幅は宇宙線加速の要だが、その種磁場の生成・維持機構は未解明。
2. 無衝突領域では Weibel/フィラメント化不安定が逆流プラズマの自由エネルギーを磁場に変換し、相対論(例: GRB)・非相対論の両衝撃で主要な磁場生成機構として確立 [Weibel 1959; Medvedev & Loeb 1999]。
3. SNR 衝撃では反射イオンが上流プラズマに逆流 → ビームフィラメント化で前駆体が準周期的な電流フィラメントに組織化。
4. これらの成長・併合・磁場構造形成は、フィラメントの線形安定性（主に FMI と drift-kink 枝）が支配。
5. しかし**自己無撞着な**SNR フィラメント平衡（古典 drift-kink 研究の理想 Harris シートでなく）の線形スペクトルは未特性化 ← 本研究の gap。

## 論点の取捨（今回の設計）
- **入れた**: WI が磁場生成の本命であること、SNR での逆流＝反射イオン、フィラメント→FMI/DKI、そして「非 Harris の自己無撞着平衡の線形解析が未整備」という gap。
- **省いた**: WI の歴史的詳細（anisotropy-driven vs beam-driven の区別）、宇宙線加速理論の詳細、磁化・ガイド磁場の話。abstract の背景としては上記5点で十分で、詰め込みすぎない。
- gap の一文（案A第5文）が、後段の結果（Pritchett 非 MHD drift kink との一致・Ng2019 kink 判定）へ自然につながる。

## 引用リスト（2026-07-28 検証結果）
検証済（~/paper で全文/要旨確認 ★）:
- **Weibel, E. S. (1959)** *Phys. Rev. Lett.* 2, 83. — Weibel 不安定の原典（well-known, 手持ちなしだが定番）。
- **Medvedev, M. V. & Loeb, A. (1999)** *ApJ* 526, 697. — 相対論二流体(Weibel)が GRB 無衝突衝撃で磁場生成 ★。
- **Silva, L. O. et al. (2003)** *ApJ* 596, L121. — 衝突プラズマ殻の PIC、near-equipartition Weibel 磁場生成。応用先に SNR を明示 ★。
- **Matsumoto, Y., Amano, T., Kato, T. N. & Hoshino, M. (2015)** *Science* 347, 974（doi:10.1126/science.aaa3145）。「Stochastic electron acceleration during spontaneous turbulent reconnection in a strong shock wave」。**supercritical collisionless shock でイオンが部分反射され上流で counterstream し kinetic 不安定を励起**＝我々の反射イオン→フィラメント像そのもの ★。**第3文の主引用に最適**（タイトルページ画像で確認; pdftotext は綴じ込みの隣接化学論文を誤抽出していた）。
- **Kato, T. N. (2007)** *ApJ* 668, 974. — 無磁化 e⁻e⁺ の Weibel 衝撃。下流粒子が上流へ戻り上流で電流フィラメント生成 ★。ただし**ペアプラズマ・相対論**なので補助的に。

要確認・注意:
- 非相対論イオンの Weibel 衝撃をさらに強調するなら Kato & Takabe (2008) *ApJ* 681, L93 等も候補（手持ちなし、要確認）。

後段（結果）で引く drift-kink / 線形理論（別途確認済）:
- Pritchett, Coroniti & Decyk (1996) *JGR* 101, 27413（要旨のパリティ命題★）; Daughton (1998) *JGR*; Ng, Hakim, Juno & Bhattacharjee (2019) *JGR* 124（Fig.2 ★）。

## 推奨引用配置（背景パート）
- 第2文（磁場生成）: `[Weibel 1959; Medvedev & Loeb 1999]`（＋任意で Silva 2003）。
- 第3文（反射イオン→上流フィラメント）: `[Matsumoto et al. 2015]`（主）＋任意で `[Kato 2007]`。Matsumoto+2015 が supercritical shock の反射イオン counterstream を明示するので最適。

---

## 確定版 v3（2026-08-07, SGEPSS 予稿 — 教授指摘を反映）

> 教授指摘（2026-08-07）を反映した全面改稿版。主な変更: ①全体に参考文献を挿入、②導入で「磁場構造の解明が粒子加熱・加速と磁場観測の説明に効く」動機を明示、③「電流層ができる→壊れる（リコネクション/キンク）→崩壊＝粒子加熱・散逸」の流れを追加、④Vanthieghem+2018（相対論的FMI/DKI競合）を挙げ「その非相対論・SNR版が本研究」と位置づけ、⑤できたことを2段階（平衡構築／線形解析）に分けて明示、⑥今後の課題（パラメータ依存・他の波数方向・非線形シミュレーション）を最終段落に。細部（成長率の数値、包絡線の分岐構造など）は書かない方針。

### 引用文献 — 引用意図と根拠箇所（2026-08-07, 全て ~/paper/ の PDF 本文で確認）

**Medvedev, M. V. & Loeb, A. (1999)** *ApJ* **526**, 697.
"Generation of Magnetic Fields in the Relativistic Shock of Gamma-Ray Burst Sources"
- 引用意図: WI（相対論的二流体不安定）が無衝突衝撃波で磁場を生成する機構であることの基礎文献。
- 根拠箇所: Abstract（p.697）「We show that the relativistic two-stream instability can naturally generate strong magnetic fields with 10⁻⁵–10⁻¹ of the equipartition energy density, in the collisionless shocks of gamma-ray burst (GRB) sources.」
- ⚠️ **注意**: 本論文の対象は **GRB の相対論的衝撃波**であり、SNR でも非相対論的高マッハ数でもない。現行の本文はこの文献を「SNR衝撃波のようにマッハ数が非常に高い場合には」という文に付けており、**対象が一致していない**。下記「未対応事項」参照。

**Matsumoto, Y., Amano, T., Kato, T. N. & Hoshino, M. (2015)** *Science* **347**, 974.
"Stochastic electron acceleration during spontaneous turbulent reconnection in a strong shock wave"
- 引用意図: ①反射イオンと入射イオンの逆流を自由エネルギーとする ion beam Weibel が電流層を作ること（＝本研究の3粒子種平衡そのものの描像）、②その電流層上で自発的な磁気リコネクションが起きること、③リコネクションがエネルギー散逸・粒子加速の担い手であること。
- 根拠箇所:
  - p.977 本文「Figure 2 demonstrates how the current sheets form spontaneously. In the supercritical shock, a portion of the upstream ions (∼20%) is reflected specularly by the shock front and gyrates in the transition region … The large dispersion both in the x and z (out-of-plane) velocity components between the incident and reflected ions is the free energy for exciting the ion beam Weibel instability」
  - Fig. 2 キャプション「Mechanism of current sheet formation … Illustration of the current sheet formation mechanism via the ion beam Weibel instability」
  - Abstract「efficient electron energization can occur during turbulent magnetic reconnection arising from a strong collisionless shock」「magnetic reconnection as an agent of energy dissipation and particle acceleration in strong shock waves」
  - p.977「the Alfvén Mach numbers are usually much higher than 10 in supernova remnant shocks」← SNR の高マッハ数を述べているのはこちら

**Jikei, T., Amano, T. & Matsumoto, Y. (2024)** *ApJ* **961**, 157.
"Enhanced Magnetic Field Amplification by Ion-beam Weibel Instability in Weakly Magnetized Astrophysical Shocks"
- 引用意図: ①WI が作った電流層上で自発的リコネクションが起きること、②それが若い SNR 衝撃波に適用できること、③フィラメントがキンク不安定性で崩壊しうること。
- 根拠箇所:
  - Abstract「Particle-in-cell simulations for magnetized electrons identify a dynamo-like mechanism of magnetic field amplification, which eventually leads to spontaneous magnetic reconnection. We conclude that this scenario is applicable to typical young supernova remnant shocks.」
  - §1 Introduction「the current filaments might eventually break up by the kink instability in a realistic fully 3D system」
  - §5 Discussion「the unmagnetized Weibel instability, in which large-scale filamentary currents are generated and eventually disrupted via the kink instability (Ruyer & Fiuza 2018)」

**Pritchett, P. L., Coroniti, F. V. & Decyk, V. K. (1996)** *JGR* **101**(A12), 27413.
"Three-dimensional stability of thin quasi-neutral current sheets"
- 引用意図: ①薄い電流層の崩壊機構として drift-kink モードが知られていること、②本研究のパリティ判定基準の出典、③tearing も同じ系で不安定化すること。
- 根拠箇所: Abstract（p.27413）
  - 「the drift kink mode is found to be of critical importance」
  - 「the drift kink mode is a non-MHD mode with a polarization structure such that E1y is an antisymmetric function of z while E1z is a symmetric function with E1z(0) ≠ 0」← **本研究の「軸方向電場の中性線に対するパリティ」判定の根拠**
  - 「The kyL∼1 drift kink modes are always the first to grow in the simulations; subsequently, tearing-like modes with a dominant kx wave vector also become unstable.」← **tearing の引用根拠にもなる**

**Vanthieghem, A., Lemoine, M. & Gremillet, L. (2018)** *Phys. Plasmas* **25**, 072115.
"Stability analysis of a periodic system of relativistic current filaments"
- 引用意図: ①filament merging (FMI) が電流フィラメント系の崩壊機構として知られていること、②相対論的ペアプラズマで FMI と DKI の競合を扱った先行研究であり、本研究はその非相対論・イオン電子版という位置づけ。
- 根拠箇所:
  - Abstract「we investigate the stability of a stationary periodic chain of nonlinear current filaments in counterstreaming pair plasmas … a weakly nonlinear symmetric system, prone to purely transverse merging modes; a strongly nonlinear symmetric system, dominated by coherent drift-kink modes … we derive an analytical criterion for the transition between the dominant filament merging and drift-kink instabilities in symmetric two-beam systems」
  - §VI 結論「as the nonlinearity of the equilibrium filaments … increases, our numerical Floquet-type calculations predict a smooth transition from a dominant, purely transverse filament merging instability (FMI) to a comparatively long-wavelength drift-kink instability (DKI) … FMI is progressively mitigated and finally overtaken by DKI modes when ξ ≳ 2.5」

### 日本語

　無衝突衝撃波は長年にわたり研究されてきた。その主な動機は、粒子加速機構の理解である。増幅された磁場は粒子を散乱し、加熱と加速の効率を決める。また、観測から推定される磁場強度を説明する機構としても注目されてきた。したがって、この磁場がどのように生成され、どのような構造をとるかを明らかにすることが本質的である。無衝突衝撃波では、Weibel instability (WI) が磁場を生成する機構となることが示されている [Medvedev & Loeb, 1999]。とくに超新星残骸 (SNR) 衝撃波では Alfvén マッハ数が 10 を大きく超え [Matsumoto et al., 2015]、WI が磁場生成の支配的な機構と考えられている。

　WIが衝撃波遷移層に電流フィラメント構造を形成することは、すでに確立している。重要なのはその後の発展である。フィラメント間に形成された電流層自体が不安定となり、崩壊する。この崩壊によって、磁場のエネルギーが粒子の加熱として散逸する。SNR衝撃波を模したPICシミュレーションでは、これらの電流層上で自発的な磁気リコネクションが生じることが示されており [Matsumoto et al., 2015; Jikei et al., 2024]、キンク型の崩壊も示唆されている。崩壊を引き起こす不安定性の候補としては、tearing、drift-kink [Pritchett et al., 1996]、filament merging [Vanthieghem et al., 2018] が知られている。

　Vanthieghem et al. [2018] は、相対論的ペアプラズマにおける filament merging と drift-kink の競合を調べた。本研究では、これをSNR衝撃波に対応する非相対論的なイオン・電子系で行うことを目指す。この種の電流層に対する既存の線形解析は、実際に形成される自己無撞着な状態ではなく、Harris sheet のような理想化された平衡に依拠している。そのため、これらのモードのうちどれが実際に生じるのかは明らかでない。

　本研究は2段階で進める。第一に、電子・入射イオン・反射イオンの3流体モデルを用いて、周期的な平衡解を数値的に構築する。WIがある程度成長してフィラメント化した状態を、無摂動の平衡状態とみなす。第二に、この平衡のまわりで線形固有値問題を解き、tearing、drift-kink、filament merging、sausage のうちどのモードが現れるかを調べる。

　代表的な結果として、電流に平行な波数ベクトルに対し、kinkモードとsausageモードの双方が得られた。これらは、軸方向電場の中性線に対するパリティによって区別される。

　今後の課題は複数ある。まず、線形スペクトルのパラメータ依存性、とくにフィラメントのピンチの強さに応じてどのモードが優勢になるかは、まだ明らかにしていない。また、斜め方向に伝播するモードが支配的となる可能性があるため、電流に平行でない波数ベクトルへの拡張も必要である。ただし、最終的にどのモードが発達して電流層を崩壊させるかは、非線形段階でしか決まらない。この点については、同じ平衡を初期条件とするPICシミュレーションによって調べる予定である。

### English

　　Collisionless shocks have been studied for decades, largely motivated by the goal of understanding particle acceleration. The amplified magnetic field scatters particles and sets the efficiency of their heating and acceleration, and has also drawn attention as the mechanism accounting for the field strengths inferred from observations. Clarifying how this field is generated, and what structure it takes, is therefore essential. At collisionless shocks, the Weibel instability (WI) has been shown to generate magnetic field [Medvedev & Loeb, 1999]. At supernova remnant (SNR) shocks in particular, where the Alfvén Mach number greatly exceeds 10 [Matsumoto et al., 2015], WI is regarded as the dominant field-generation mechanism.

　　That WI forms a current-filament structure in the shock transition layer is by now established. What matters is the subsequent evolution. The current sheets formed between these filaments become unstable and break up, and through this breakup the magnetic energy is dissipated as particle heating. Particle-in-cell simulations modeling SNR shocks have shown spontaneous magnetic reconnection occurring on these sheets [Matsumoto et al., 2015; Jikei et al., 2024], and kink-type disruption has also been suggested. Known candidates for the instability driving this breakup are the tearing, drift-kink [Pritchett et al., 1996], and filament-merging [Vanthieghem et al., 2018] modes.

　　Vanthieghem et al. [2018] examined the competition between filament merging and drift kink in a relativistic pair plasma. Here we aim to carry this out for the non-relativistic, ion-electron system corresponding to SNR shocks. Existing linear analyses of such current sheets rely on idealized equilibria such as the Harris sheet, rather than the self-consistent state that actually forms. It therefore remains unclear which of these modes actually arises.

　　We proceed in two stages. First, using a three-fluid model — electrons, incoming ions, and reflected ions — we construct a periodic equilibrium numerically, regarding the filamented state reached after WI has grown as the unperturbed equilibrium. Second, we solve the linear eigenvalue problem about this equilibrium to determine which of the tearing, drift-kink, filament-merging, and sausage modes appear.

　　As a representative result, for wavevectors parallel to the current we obtained both kink and sausage modes. These are distinguished by the parity of their axial electric field about the neutral lines.

　　Several directions remain. First, the parameter dependence of the linear spectrum — in particular which mode becomes dominant as the filament pinching is varied — has yet to be established. The analysis also needs to be extended to wavevectors not parallel to the current, since obliquely propagating modes may become dominant. Which mode ultimately develops and disrupts the current sheets, however, can only be settled in the nonlinear stage; we plan to address this with particle-in-cell simulations initialized from the same equilibrium.

### v3 の未対応事項
- **英語本文は 2733 字（空白込, 403 語）で、2000 字制限を 733 字超過**。提出前に要削減。段落別の内訳は para1=572 / para2=618 / para3=427 / para4=405 / para5=184 / para6=527。削減候補は para2（リコネクション文献の列挙を圧縮）と para3（Vanthieghem の説明を1文に）。
- ✅ **Medvedev & Loeb (1999) の引用位置の不一致は 2026-08-07 に修正済**。同論文は GRB の相対論的衝撃波が対象で SNR を扱っていないため、文を2つに分割した: 「無衝突衝撃波では WI が磁場を生成する機構となることが示されている [Medvedev & Loeb, 1999]。とくに SNR 衝撃波では Alfvén マッハ数が 10 を大きく超え [Matsumoto et al., 2015]、WI が支配的な機構と考えられている。」SNR の高マッハ数は Matsumoto+2015 p.977 が根拠。
- tearing 単独の一次文献は引いていないが、Pritchett+1996 の Abstract に「subsequently, tearing-like modes with a dominant kx wave vector also become unstable」とあり、drift-kink と tearing の両方の根拠になる。追加引用（Furth+1963 等）は必須ではない。
- 「代表的な結果」の段落は、包絡線が複数ブランチの乗り換えを含むこと（[[snr-current-filament-linear-analysis]] 関連、2026-08-07 のブランチ追跡で判明）には触れていない。細部を書かない方針のため意図的に省略。

---

## 旧版 v2（2026-08-07, 教授指摘反映前）

> 案A/B（電子加速→SNR）とは動機の立て方が異なる版（電子加速問題→イオンビーム不安定性→WI→高マッハ数のSNR、という流れ）。v1（背景+導入のみ）に、研究紹介パート（平衡解の構成・妥当性確認の方針・これまでの結果・今後）を追加した全体版。引用は未挿入（保留中）。教授承認済みなのは warm natural 平衡の線形解析結果のみ（[[snr-abstract-scope-natural-only]]）、cold平衡の数値は書かないこと。PICとの成長率比較は warm natural については未検証なので本文に書いていない。

### 日本語

　無衝突衝撃波は長年にわたり研究されてきた。その主な動機は、荷電粒子の加速機構の理解である。中でも、衝撃波遷移層におけるイオンビーム不安定性は、電子加速問題を解く鍵として最も重要視されてきた。超新星残骸 (SNR) のような高マッハ数領域では、Weibel instability (WI) が支配的なモードとして注目されている。

　SNRの衝撃波領域では、WIに関するシミュレーション研究が行われてきた。それらの研究では、WIによって形成された電流層上でリコネクションなど他のモードが成長することが示唆されている。しかし、既存の電流層の線形解析はハリスシートのような理想的な平衡状態に依拠しており、この電流層の平衡状態に基づいた線形解析は行われていない。したがって、シミュレーションで観測されるこれらのモードが線形不安定性と一致するか、また他にどのような不安定モードが存在しうるかは明らかになっていない。

　本研究では、SNR衝撃波遷移層に対して周期的な電流フィラメント構造の平衡解を構築し、その線形安定性解析によって、どのようなモードが存在するかを解明する。

　この平衡解は電子・反射イオン・入射イオンの3粒子種で構成され、SNR衝撃波遷移層で想定される電流フィラメント構造を模している。ただし、実際のSNRのパラメータでは平衡解の形が歪む。そこで解析の妥当性を確立するため、まず平衡解が滑らかなままとなるパラメータ領域で固有値問題を解く。

　波数ベクトルが電流と平行な摂動については、最速成長モードは軸方向電場の中性線に対するパリティ（奇）で同定されるkink型不安定性である。軸方向波数を系統的に変えながら調べた結果、識別したkink・sausage以外にこれらより速く成長するモードは存在しないことを確認した。現在、歪んだSNR実パラメータの平衡、およびフィラメントが孤立Harris sheetに近づく強ピンチ極限へと解析を拡張し、ここで見つかったモード構造を古典的なdrift-kink不安定性に接続することを進めている。

### English

　　Collisionless shocks have been studied for decades, largely motivated by the goal of understanding the acceleration of charged particles. In particular, Ion-beam instabilities in the shock transition layer have long been considered central to the electron acceleration problem. At a very high Mach number regime, such as in supernova remnant (SNR), the Weibel instability (WI) has drawn attention as the dominant mode.

　　Simulations of WI have been carried out for SNR shocks. These studies suggest that other modes, such as reconnection, grow on the current sheets formed by WI. However, existing linear theory for such current sheets relies on idealized equilibria, such as the Harris sheet, rather than the self-consistent state that actually forms in this system. It therefore remains unclear whether these modes correspond to linear instabilities of the equilibrium, and what other modes may exist.

　　In this study, we construct an equilibrium solution for a periodic current-filament structure in the SNR shock transition layer. We then perform its linear stability analysis to determine which modes are present.

　　This equilibrium comprises three particle species — electrons, reflected ions, and incoming ions — and models the current-filament structure expected in the SNR shock transition region. With SNR-realistic parameters, however, the equilibrium profile becomes distorted. To establish the validity of our analysis, we therefore first solve the eigenvalue problem in a parameter regime where the equilibrium remains smooth.

　　For perturbations whose wavevector is aligned with the current, the fastest-growing mode is a kink-type instability, identified by the odd parity of its axial electric field about the neutral lines. A systematic scan over the axial wavenumber confirms that no branch grows faster than the kink and sausage modes we identify. We are now extending the analysis to the distorted, SNR-relevant equilibrium, and tracking how the growth rate and eigenfunction of this mode change as the filaments are pinched toward isolated Harris sheets, where classical drift-kink predictions are available for direct comparison.

### 未対応の指摘（次回持ち越し）
- "Ion-beam" の大文字I（文中なので小文字に）。
- 英語第4段落の "shock transition region" と、第1・3段落の "shock transition layer" が不統一（layer に揃えるのを推奨）。
- 英語第5段落の "For perturbations whose wavevector is aligned with the current" は冗長（"For wavevectors parallel to the current," への短縮を提案済み、未反映）。
- 第1段落の「電子加速問題を解く鍵」という動機の強さと、第3段落の結論（モード同定どまり）の接続をどうするか（弱める案を提示済み、保留中）。
- 引用は未挿入。提出前に上記「推奨引用配置」を参考に挿入。
- warm natural 平衡の線形解析結果とPIC成長率の比較は未実施。実施しない限り、比較に言及する文を本文に加えないこと。
