# DKI/FMI 固有関数プロファイル：文献（Harris シート）との比較と前提差

教授MTG(7/27)指針B：DKI/FMI 固有関数を既存論文（多くは Harris シート平衡）と比較し、前提条件の違いを考慮する。
本ノートは (1) 文献の DKI/FMI 固有関数の典型構造、(2) 文献 vs 我々の平衡の前提差、(3) 我々の warm natural モードの位置づけ、をまとめる。

> **引用の扱い（citation-verification）**：以下の書誌は Web 検索で publisher/arXiv ページを確認したもの。全文を開けたのは Gingell+2014（arXiv 全文）のみ。他は書誌情報とアブストラクト要約に基づく（★=全文確認、☆=書誌+要約のみ）。abstract に載せる前に該当箇所を原著で最終確認すること。

---

## 1. 文献の DKI（ドリフト・キンク）固有関数の典型構造

### 平衡（ほぼ全て Harris シート）
古典的 DKI 研究は**単一 Harris 電流シート**を平衡にする（例、Gingell+2014 ★ の式）：

```
B_y(x) = B0 sin(θ/2) tanh(x/L)      （反平行 θ=π で B0 tanh(x/L)）
n(x)   = n0 + n_cs sech²(x/L)
```

- **ゼロ電場条件（圧力釣り合い）** `P⊥,e + P⊥,p + B²/8π = const` → `β⊥,e + β⊥,p = 1`。これは Harris 平衡の帰結で、**静電場 E=0**（Φ=0）。教授が以前指摘した「Harris では V_i/T_i+V_e/T_e 型の制約が要る」に対応。
- 電流は**片方のドリフトイオン**（または電子）が担う：`v_{p,0} = -w⊥ρ_p/L`。

### DKI モードの性質
- **k ∥ J**（電流＝ビーム方向に波数を持つ）で、シートを**横に折り曲げる（kink）／絞る（pinch=sausage）** 変位。Gingell+2014 ★ は「温度異方性がシートを pinch/kink する narrow-band drift-kink」と記述。
- 二流体論（Yoon, Lui & Sitnov 1998 ☆, JGR 103, 11875）：中性シートについて**対称・反対称両方の磁場摂動**を持つ成長モードが出る。成長率は実周波数と同程度で、シート外イオンジャイロ周波数の相当割合。
- 運動論（Daughton 1998 ☆ JGR 103; Daughton 1999 ☆ *Phys. Plasmas* 6, 1329「The unstable eigenmodes of a neutral sheet」）：シート中心付近で**イオン慣性長より細かい空間構造**、**有意な静電成分**を持ち real frequency を変える。現実的質量比・厚いシートでは DKI 成長率は tearing より弱い。
- **局在**：摂動は電流シート近傍（勾配領域）に局在。

### DKI 固有関数の具体形（Ng, Hakim, Juno & Bhattacharjee 2019 ★, 全文・図確認）
二流体 five/ten-moment で Harris シート（`B_x=B0 tanh(z/L)`, `n=n0 sech²(z/L)`, 平衡条件 `β_e+β_i=1`）の drift-kink / drift-sausage を直接固有値計算。**Fig.2 が固有関数の具体形**（\(E_y\) の Re 実線・Im 破線, `z/L`=−10..10）：

- **kink（DKI）モード ＝ \(E_y\) が中性シート(z=0)について奇（antisymmetric, odd）**。
- **sausage モード ＝ \(E_y\) が偶（symmetric, even）**。← 両者はパリティで截然と区別される。
- **局在**：固有関数は `z/L≈±5` 以内に集中し `±10` で 0 に減衰（シート近傍局在）。
- **長波長**：最速 kink は `k_y L≈0.5`（five-moment）〜`k_y ρ_i≈1`（ten-moment）＝**シート厚スケール**。
- **成長率**：`γ/Ω_ci≈0.05–0.22`（イオンジャイロ周波数単位）、`T_i/T_e` を下げると増、**質量比を上げると減**（運動論と整合）。

> **座標対応**：Ng+2019 の `z`（シート横断）＝我々の `y`（フィラメント横断）、彼らの `k_y`（電流/ドリフト方向）＝我々の `k_x`（ビーム方向）。よって「DKI 固有関数の中性シートについての奇パリティ」は、我々では「フィラメントの \(B_{0z}\) ゼロ交差線について \(E_x\)(彼らの \(E_y\) に相当) が奇」に対応する。

### 我々の設定に最も近い文献：ion-ion kink
- **Karimabadi, Daughton, Pritchett+ 2003 ☆**（JGR 108, doi:10.1029/2003JA010026）「Ion-ion kink instability in the magnetotail: 1. Linear theory」。**2つのイオン集団が相対ドリフト**する系のキンクで、電子電流由来の古典 DKI と区別される。我々の inc/ref 逆流イオンビームに構造的に最も近い。

## 2. FMI／フィラメント（Weibel系）の典型構造
- 逆流ビームのフィラメント化（Weibel/filamentation）：磁場を指数成長させ、プラズマを磁化。**k⊥beam** の純フィラメント（我々の kx=0）と、有限 kx の**斜め（oblique）モード**が連なる。
- Bret, Gremillet, Dieckmann ら ☆：波数平面上で**成長率の近い斜めモードが複数共存**し、低波数モードから高波数の斜めモードへ移り得る。← 我々の「kx 掃引で FMI→斜め枝が連続」と整合。

## 3. 前提差の表（文献 Harris-DKI vs 我々の平衡）

| 項目 | 文献の Harris-sheet DKI | 我々の SNR 平衡 | 帰結（固有関数への影響） |
|---|---|---|---|
| 配位 | **単一**電流シート | **周期フィラメント配列**（λ0 周期） | Bloch 波数 K が入り、包絡 e^{iKy} が対称性を崩す |
| 磁場 | B∼tanh(x/L)（1シート） | B0z∼周期的（正弦波状, 奇） | 局在でなく周期構造 |
| 密度 | n∼sech²(x/L) | inc/ref が緩やかに変調（空乏しない warm） | シート局在でなく全域に分布 |
| 電場 | **E=0**（β_e+β_p=1 の Harris 制約） | **E0y≠0**（Harris 制約を破る, Φ≠0） | 背景場力項が入り、モード構造・FMI 抑制が変わる |
| 電流の担い手 | 片方のドリフトイオン（or 電子） | **2本の逆流イオンビーム**(inc/ref) | 駆動は ion-ion 型＝Karimabadi2003 に近い |
| 駆動源 | シート電流＋（温度異方性） | ビーム間相対ドリフト（自由エネルギー） | 二流体/フィラメント的 |
| 磁化 | 有限ガイド磁場ありが多い | σ=0（非磁化ベース） | O偏波分離、P偏波のみ |
| モード | k∥J のキンク、δB はシートに局在・sym/antisym | 有限 kx の斜め枝、δB は周期内に微細構造 | パリティ・局在が本質的に異なる |

## 4. 我々の warm natural ピークの位置づけ

- 前段の固有関数診断（別図 `eigfun_convergence.png`）で、warm ピーク(kx/k0=1.10)は **DKI キンクに特有の構造（δBz∝∂_yB0z, ビーム同位相, シート局在）を満たさない**。
- 文献対照でこれは自然に説明できる：我々のモードは**逆流イオンビーム駆動の斜め枝**で、古典的な単一 Harris シートの電子 DKI ではない。最も近い文献は **ion-ion kink（Karimabadi 2003）** と **フィラメント斜めモード（Bret ら）**。
- **固有関数の非対称性**も前提差で説明：文献 Harris-DKI が中性シートについて sym/antisym の確定パリティを持てるのは、単一シート＋E=0 の対称平衡だから。我々は (i) 周期配列＋Bloch e^{iKy}、(ii) 有限 kx の伝搬、(iii) E0y≠0（Harris 制約破れ）により確定パリティを持たない。

### 具体的な照合基準（Ng+2019 Fig.2 を物差しに）
DKI と呼べるかは、次の**具体的な固有関数の特徴**で判定できる（Ng+2019 の DKI kink モードの特徴）：
1. 電流方向電場（我々の \(E_x\)）がフィラメント \(B_{0z}\) ゼロ交差線について**奇（antisymmetric）**。
2. **シート近傍に局在**（周期全体に拡散しない）。
3. **長波長**（`k_x·L`〜0.5–1 = シート厚スケール）。

我々の warm ピーク(kx/k0=1.10)はこれらを満たさない（高 kx・非局在・確定パリティなし）。→ **次の具体作業**：我々の**低 kx（シート厚スケール）**の枝を取り出し、単一フィラメントの中性線について \(E_x\) の奇パリティ・局在を Ng+2019 と同じ物差しで確認する。あれば「DKI 的枝が存在」と言え、なければ「DKI 枝は立たない（純粋に斜めビーム不安定）」と結論できる。

### 検証結果：K=0（純 k∥J）で kink 枝を同定（2026-07-27）
純 DKI は k∥J＝横方向波数ゼロ。我々でこれは **K=0**（フィラメント間の包絡波数ゼロ）に対応する。K=0 では摂動が純周期になり
パリティが確定するので、Ng+2019 と同じ物差しで kink/sausage を分離できる（駆動 `work/python/plot_kink_sausage_K0.py`）。

- **kink（E_x 奇）と sausage（E_x 偶）が K=0 で厳密分離**（混合ゼロ）。E_x 奇 ⟺ B_z 偶が一貫。
- **低 kx（kx/k0≲0.35）で kink 枝が sausage より明確に優越**（例 kx/k0=0.15 で kink γ=9.3×10⁻⁴ ≫ sausage 2.2×10⁻⁶）＝**DKI 的レジームが実在**。高 kx（≳0.45）では kink≈sausage で縮退（2本の中性線が独立化）。
- **代表 kink（kx/k0=0.15）の固有関数**（図 `kink_eigenfunction_K0.png`）は Ng+2019 の kink 判定を満たす：**δE_x が中性線 (B0z=0) について奇**、背景 B0z をなぞる、**δB_z は偶で大振幅**、滑らかな単一構造（高 kx ピークの微細構造と対照的）。
- **ただし従属**：kink（DKI）枝の γ≈9×10⁻⁴ は高 kx ピーク（kink/sausage 縮退の斜め枝, γ≈7.9×10⁻³）の**約 1/8**。つまり我々の系で DKI は「存在するが最速ではない低 kx の従属枝」。

### Pritchett 1996 の偏波命題との照合（2026-07-28）＝ 完全一致
Pritchett, Coroniti & Decyk (1996, *JGR* 101, 27413) の**要旨**は、初期 Bz なしの drift kink を
「**非 MHD モード**で、\(E_{1y}\)（電流方向）は奇（antisymmetric）、\(E_{1z}\)（横断方向）は偶（symmetric）で
\(E_{1z}(0)\neq0\)（中性面で非ゼロ）」と述べる。\(E_{1z}(0)\neq0\) が**非 MHD の署名**（理想 MHD なら
\(E_{1z}=v_{1y}B_0/c\propto B_0\to0\) で中性線で消えるが、drift kink は消えない）。
座標対応：Pritchett \(E_{1y}\)＝我々 \(E_x\)、\(E_{1z}\)＝我々 \(E_y\)。我々の K=0 kink（kx/k0=0.15）で照合（<code>work/python/plot_pritchett_parity_compare.py</code>）：

- **\(E_x\)（電流方向）: 奇（odd=1.00）・中性線で 0** → Pritchett \(E_{1y}\) と<strong>一致 ✓</strong>
- **\(E_y\)（横断方向）: 偶（even=1.00）** → Pritchett \(E_{1z}\) の対称性と<strong>一致 ✓</strong>
- **\(E_y\) は中性線で非ゼロ**（\(|E_y|/\max\)=0.82（y=0）, 1.00（λ0/2））→ Pritchett \(E_{1z}(0)\neq0\) と<strong>一致 ✓</strong>

**3点すべて一致**。我々も \(E_{0y}\neq0\)・warm で非 MHD なので、\(E_y\) が中性線で消えない＝Pritchett の
**非 MHD drift kink と同じ偏波構造**。これは前提差（周期配列・非 Harris）を越えて成り立つ頑健な一致で、
我々の低 kx・K=0 モードが drift kink であることの、Pritchett に基づく強い裏付けとなる。

### 結論（abstract 向け）
- 「DKI」と単純に名指すのは前提差が大きく危うい。**「逆流イオンビームによる有限 kx の斜め（ion-ion kink / oblique filamentation 系）不安定」**と記述するのが文献整合的で安全。
- 比較の枠組み：Harris-sheet DKI（Daughton/Yoon/Gingell/Ng+2019）は参照点として引きつつ、**我々の平衡は非 Harris（周期配列・E0y≠0・ion-ion 駆動）**である点を明記して差分を説明する。
- 具体照合は Ng+2019 Fig.2 の kink モード特徴（\(E_x\) 奇・シート局在・長波長）を物差しにする。

## 参考文献（要最終確認）
- Daughton, W. (1998) "Kinetic theory of the drift kink instability in a current sheet", *JGR* 103, doi:10.1029/1998JA900028. ☆
- Daughton, W. (1999) "The unstable eigenmodes of a neutral sheet", *Phys. Plasmas* 6, 1329. ☆
- Yoon, P. H., Lui, A. T. Y., Sitnov, M. I. (1998) "Two-fluid theory of drift-kink instability in a one-dimensional neutral sheet", *JGR* 103, 11875. ☆
- Karimabadi, H., Daughton, W., Pritchett, P. L. (2003) "Ion-ion kink instability in the magnetotail: 1. Linear theory", *JGR* 108, doi:10.1029/2003JA010026. ☆
- **Ng, J., Hakim, A., Juno, J., Bhattacharjee, A. (2019)** "Drift instabilities in thin current sheets using a two-fluid model with pressure tensor effects", *JGR Space Physics* 124, doi:10.1029/2018JA026313（arXiv:1903.09618）. **★（全文・Fig.2 確認）** ← DKI 固有関数の具体形の主典拠。
- Gingell, P. W., Burgess, D., Matteini, L. (2014) "The three-dimensional evolution of ion-scale current sheets: tearing and drift-kink instabilities...", arXiv:1411.4422. ★
- Bret, A., Gremillet, L., Dieckmann, M. E. ら：filamentation/Weibel と斜めモードの比較（要書誌確定）. ☆
