# SNR Filament Stability

SNR前駆体衝撃波を模した3成分（電子・入射イオン・反射イオン）電流フィラメント平衡の構築と、線形安定性解析。

Vanthieghem et al. (2018) のペアプラズマ平衡アプローチを非相対論的な3成分系に拡張し、drift-kink instability (DKI) と filamentation merging instability (FMI) の競合を Fourier-Bloch 法で調べる。

## できること

- **平衡構築**: 電子静止系での周期的電流フィラメント平衡を数値的に解く（`src/fmi/snr_equilibrium.py`）
- **線形安定性解析**: Fourier-Bloch 法による固有値問題 ω x = L(kx, K) x を解き、DKI/FMI の成長率と固有関数を求める（`src/fmi/snr_linear_bloch.py`）
- **PIC比較**: マッハ数ベースのパラメータ設計（`src/fmi/mach_parameters.py`）と Weibel 加熱後温度（`src/fmi/wi_saturated_temperature.py`）で、PIC シミュレーションと直接比較できる平衡を用意する

## ディレクトリ構成

- `src/fmi/`: 平衡ソルバー・線形演算子・パラメータ設計の本体コード
- `notebooks/`: marimo インタラクティブノートブック
- `tests/`: 自動テスト
- `docs/`: 解析結果と導出のHTML/Markdown資料（GitHub Pagesで公開）
- `work/`: 図・CSVなどの生成物（gitで追跡しない）

## クイックスタート

```bash
uv sync
uv run pytest
```

## 公開ページ

- Stable (`main`): <https://ishida-kohki.github.io/snr-filament-stability/main/>
- Development (`develop`): <https://ishida-kohki.github.io/snr-filament-stability/develop/>

## ブランチ運用

- `main`: 安定版
- `develop`: 開発中の作業（通常はここにコミット）
- `feature/*`: 個別の実験・機能追加
