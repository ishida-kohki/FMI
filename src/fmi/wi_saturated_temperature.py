"""Weibel 不安定性で加熱された準平衡状態の温度（PIC 実測値）。

背景
----
`mach_parameters.build_solver_params` が返す温度は**上流の**一般的な SNR 環境値
（Te=0.010, Tinc=Tref=6.25e-04）である。しかし本研究が 0 次解として扱うのは
「WI がある程度成長した結果として生じる準平衡状態」なので、その時点の温度は
上流値より上昇している。上流値をそのまま使うと、平衡ソルバの Boltzmann 因子
n_s ∝ exp(-q_s(Φ-β_s A)/T_s) が a0=0.25 で n_ref = 6e-54 という非物理的な値を返し、
線形解析が破綻する（docs/density_regularization.md §4.5）。

出典
----
`~/pic-nix-weibel/pic/example/foot/weibel/data_sigma0_m400`
（σ=0 非磁化コントロールラン、mime=400、Nx=8, Ny=512, alpha=0.2、完全版）を
`temperature_verification.py` の熱圧力定義で解析した実測値。ω_pi t=70 は
初回飽和窓（eqcmp と同じ比較点）。t=0 では PIC 温度が入力値の 0.94-0.99 倍で
一致しており、単位変換（自種質量規格化 → m_e c² 規格化、×[1, mime, mime]）の
正しさも確認済み。

加熱は強い異方性を持つ（Weibel 変調方向 y が卓越）。反射イオンで T_y/T_x = 8.4。
平衡の Boltzmann 因子は y 方向の閉じ込めから来るので、等方版ではなく **T_y 版**が
物理的に対応する。実際、T_y 版でのみ線形解析が解像度収束し、正則化非依存の
磁気モードが得られる（等方版 a0=0.45 は非収束）。

準平衡窓の定義と検証
--------------------
「準平衡」は状態の性質ではなく γ との相対的な性質である。線形解析が成立する
条件は温度が止まっていることではなく、求めたいモードが e-folding する間に
背景が動かないこと:

    |d ln X/dt| < γ/10     （X = δB_z rms, λ, 各種 T_iso と T_y）

変化率は窓幅 1/γ での ln X 最小二乗傾きで測る（work/python/
quasi_equilibrium_window.py、図 work/quasi_equilibrium_window.png）。

γ = 9.1e-3 ω_pe、mime=400 より ω_pe/ω_pi=20 なので γ = 0.182 ω_pi、
閾値は 0.0182 ω_pi。判定結果:

    準平衡窓 = ω_pi t 65.5 .. 80.0（ラン終端まで、長さ 14.5 = 2.6/γ）
    ω_pi t=70 は窓の内側。最悪変化率 0.0070 ω_pi = γ/25.9。
    律速は反射イオンの T_y。λ は ω_pi t≈43 以降一定（窓内で merging なし）。

したがって従来から使っていた ω_pi t=70 はこの基準を満たす。電子と入射イオンの
温度が単調上昇中でも、上昇が γ に対して十分遅いので問題ない。

循環論法の確認: 上の γ 自体がこの t=70 の温度で得た値である。ただし判定が
覆るのは γ < 10 × 0.0070 ω_pi = 3.5e-3 ω_pe の場合に限られ、PIC 実測帯
6.6e-3〜9.0e-3 ω_pe とは 2 倍近い余裕がある。

限界: 窓はラン終端 ω_pi t=80 まで続いており、その先で閉じるかは確認できない。

正規化は snr_equilibrium と同一（k_B T/(m_e c²)）。
"""

from __future__ import annotations

__all__ = [
    "UPSTREAM",
    "WI_SATURATED_ISOTROPIC",
    "WI_SATURATED_TY",
    "PIC_SOURCE",
    "SNR_DRIFTS",
]

PIC_SOURCE = (
    "~/pic-nix-weibel/pic/example/foot/weibel/data_sigma0_m400 "
    "(sigma=0, mime=400, Nx=8, Ny=512, alpha=0.2), omega_pi*t=70"
)

# 上流（一般的な SNR 環境）: build_solver_params(MachConfig(M_S=154.92,
# beta_sh=0.25, beta_te=0.1, sigma=0, eta=0.2, mime=400)) と一致する。
# PIC の初期条件 vte=0.1, vti=vtr=0.00125 とも一致。
UPSTREAM = {"Te": 1.0000e-02, "Tinc": 6.2500e-04, "Tref": 6.2500e-04}

# WI 飽和後の等方温度 T_iso = (T_xx+T_yy+T_zz)/3。上流比 1.25 / 1.40 / 5.08。
WI_SATURATED_ISOTROPIC = {"Te": 1.2465e-02, "Tinc": 8.7457e-04, "Tref": 3.1721e-03}

# WI 飽和後の T_yy（Weibel 変調方向）。上流比 1.45 / 2.09 / 12.77。
# 平衡の y 構造に対応する温度はこちら。
WI_SATURATED_TY = {"Te": 1.4507e-02, "Tinc": 1.3089e-03, "Tref": 7.9791e-03}

# 対応するドリフト（build_solver_params の出力、自然系）。温度と組で使う。
SNR_DRIFTS = {
    "eta": 0.2,
    "beta_inc": -0.05047318611987382,
    "beta_ref": 0.19814241486068115,
}
