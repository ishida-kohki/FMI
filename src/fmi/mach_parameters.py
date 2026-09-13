"""マッハ数を代表パラメータとした前向きパラメータ設計器。

代表パラメータ (M_A, M_S) からソルバー入力と PIC config を一意に導出する。
正規化は snr_equilibrium と同一（長さ c/ω_pe、温度 T̃=k_B T/(m_e c²)、m_e=1, m_i=mime）。

速度記法（c=1、全速度を c で正規化）:
    すべての速度は c で割った無次元量とし、記号 β_<名> で表す（β_<名> ≡ <速度>/c）。
        β_sh   衝撃波速度 V_sh/c          β_A,i イオンAlfvén速度 v_A,i/c
        β_inc/β_ref/β_e  各成分ドリフト    β_th,s 熱速度 v_th,s/c
        β_cs,s 音速 c_s/c = √γ·β_th,s
    注意: 「プラズマβ」は別物（速度比ではなく圧力比 β_plasma = 2 n k_B T /(B²/2μ0)）。
    混同を避けるため本モジュールでは圧力比のみ beta_e_plasma / betae / betai と綴る。

物理チェーン（磁化時）:
    β_A,i  = √(σ/mime)                                # イオンAlfvén速度/c
    β_sh   = M_A · β_A,i                              # 衝撃波速度/c
    β_d    = β_sh/2                                   # 相対ドリフト 2·β_d = β_sh
    β_inc  = -2η·β_d / (1-(1-2η)β_d²)                 # 鏡面反射（main.cpp と同式）
    β_ref  = +2(1-η)·β_d / (1+(1-2η)β_d²)
    β_cs,i = β_sh / M_S                               # イオン音速/c
    T̃_inc = mime · β_cs,i² / γ  ( = mime·β_th,i² )    # β_cs = √γ·β_th
    β_plasma,i = 2 β_th,i²/β_A,i² = (2/γ)(M_A/M_S)²   # イオンプラズマβ（圧力比、速度ではない）

本プロジェクトの既定は非磁化（sigma=0）:
    SNR 電流フィラメント平衡・線形解析は非磁化で構築する。Jikei+ 2024 の外部磁場は
    イオンビーム（ドリフト x）に垂直で、ドリフトに直交する平均場は成分ごとに異なる
    motional 力を生み 3成分差動ドリフトの静的周期平衡と非両立のため、平衡に外部磁場を
    入れない。非磁化では β_A が無く M_A は使えないので、速度は β_sh を直接、温度は
    音速マッハ M_S、電子温度は β_te で与える。β・T は磁化版と機械精度で一致する
    （sigma は平衡 ODE に現れないため）。磁化効果が要る段は線形演算子に面外ガイド場
    √σ を加える補正で扱う（平衡は非磁化のまま）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class MachConfig:
    """マッハ数基準のパラメータ点。

    Attributes
    ----------
    M_A: イオンAlfvénマッハ数 V_sh/v_A,i（磁化時に速度をアンカー）。非磁化(sigma=0)では
         使わないため省略可（None）。
    M_S: イオン音速マッハ数 V_sh/c_s,i（温度をアンカー）。
    eta: 反射イオン割合（= PIC の alpha）。
    mime: 質量比 m_i/m_e。
    sigma: 磁化 (v_A,e/c)² = (Ω_ce/ω_pe)²。0 で非磁化。
    gamma: 音速に使う断熱指数（既定 5/3、1D断熱なら 3 も可）。
    beta_sh: 衝撃波速度 β_sh = V_sh/c。非磁化(sigma=0)時に必須の速度アンカー。
             磁化時は無視（M_A·β_A,i から導出）。
    beta_e_plasma: 電子プラズマβ（圧力比、速度ではない。磁化時の電子温度ノブ）。
    beta_te: 電子熱速度 β_th,e = v_th,e/c を直接指定（非磁化向け代替ノブ）。
    Te_over_Tinc: 電子/入射イオン温度比による代替ノブ。
    Tref_over_Tinc: 反射/入射イオン温度比（現行ラン betar=betai → 1.0）。
    """

    M_S: float
    M_A: float | None = None
    eta: float = 0.2
    mime: float = 400.0
    sigma: float = 0.0025
    gamma: float = 5.0 / 3.0
    beta_sh: float | None = None
    beta_e_plasma: float | None = None
    beta_te: float | None = None
    Te_over_Tinc: float | None = None
    Tref_over_Tinc: float = 1.0


def _beta_shock(cfg: MachConfig) -> float:
    """衝撃波速度 β_sh = V_sh/c を返す。磁化時は M_A·β_A,i、非磁化時は cfg.beta_sh。"""
    if cfg.sigma > 0.0:
        if cfg.M_A is None:
            raise ValueError("sigma>0（磁化）では M_A の指定が必須です。")
        beta_ai = math.sqrt(cfg.sigma / cfg.mime)  # β_A,i = v_A,i/c
        return cfg.M_A * beta_ai
    if cfg.beta_sh is None:
        raise ValueError("sigma=0（非磁化）では beta_sh の指定が必須です。")
    return cfg.beta_sh


def _specular_drifts(beta_d: float, eta: float) -> tuple[float, float]:
    """鏡面反射規約のシム系ドリフト (β_inc, β_ref)。main.cpp:61-62 と同式。

    beta_d は相対ドリフトの半分 β_d = β_sh/2（2·β_d = β_sh）。
    """
    s = 1.0 - 2.0 * eta
    beta_inc = -2.0 * eta * beta_d / (1.0 - s * beta_d * beta_d)
    beta_ref = +2.0 * (1.0 - eta) * beta_d / (1.0 + s * beta_d * beta_d)
    return beta_inc, beta_ref


def _electron_temperature(cfg: MachConfig, T_inc: float) -> float:
    """電子温度 T̃_e [m_e c² 単位] を独立ノブから決める（排他に1つ指定）。"""
    knobs = [cfg.beta_e_plasma, cfg.beta_te, cfg.Te_over_Tinc]
    if sum(k is not None for k in knobs) != 1:
        raise ValueError(
            "電子温度ノブは beta_e_plasma / beta_te / Te_over_Tinc の"
            "いずれか1つだけを指定してください。"
        )
    if cfg.beta_e_plasma is not None:
        if cfg.sigma <= 0.0:
            raise ValueError("beta_e_plasma は磁化時(sigma>0)のみ有効です。")
        # β_th,e² = (v_A,e/c)²·(β_plasma,e/2) = sigma·β_plasma,e/2、T̃_e = β_th,e²（m_e=1）
        return cfg.sigma * cfg.beta_e_plasma / 2.0
    if cfg.beta_te is not None:
        return cfg.beta_te ** 2
    return cfg.Te_over_Tinc * T_inc


def build_solver_params(cfg: MachConfig) -> dict:
    """MachConfig からソルバー入力 dict を返す。

    Parameters
    ----------
    cfg: マッハ数基準のパラメータ点。

    Returns
    -------
    dict:
        solve_snr_equilibrium にそのまま渡せる
        {eta, beta_inc, beta_ref, Te, Tinc, Tref}。
    """
    beta_sh = _beta_shock(cfg)  # β_sh = V_sh/c
    beta_d = beta_sh / 2.0  # 相対ドリフトの半分 2·β_d = β_sh
    beta_inc, beta_ref = _specular_drifts(beta_d, cfg.eta)

    # イオン温度: β_cs,i = β_sh / M_S、β_th,i² = β_cs,i²/γ、T̃_inc = mime·β_th,i²
    beta_cs_i = beta_sh / cfg.M_S
    T_inc = cfg.mime * (beta_cs_i ** 2) / cfg.gamma
    T_ref = cfg.Tref_over_Tinc * T_inc
    T_e = _electron_temperature(cfg, T_inc)

    return dict(
        eta=cfg.eta,
        beta_inc=float(beta_inc),
        beta_ref=float(beta_ref),
        Te=float(T_e),
        Tinc=float(T_inc),
        Tref=float(T_ref),
    )


def _boost_velocity(beta: float, beta_frame: float) -> float:
    """速度 beta を beta_frame で動く系へ移したときの値（相対論的速度合成則）。

    β' = (β - β_frame) / (1 - β·β_frame)。非相対論極限では単純な差 β - β_frame。
    """
    return (beta - beta_frame) / (1.0 - beta * beta_frame)


def to_electron_rest_frame(
    params: dict[str, float], relativistic: bool = False
) -> dict[str, float]:
    """ソルバー入力を電子静止系へ正規化する（線形解析の出発点用）。

    電流中性条件で決まる電子ドリフト β_e = (1-η)β_inc + η β_ref を 0 にし、
    イオンドリフトを同じ系へブーストする。電子感受率がドリフトなしの素の形になり、
    分散関係から Doppler 項 (ω - k·β_e c) が消えるため線形解析が簡単になる。

    既定は非相対論ブースト（単純差 β - β_e）。平衡ソルバー自体が非相対論的
    Boltzmann 近似で、β_e を線形関係 (1-η)β_inc+η β_ref から再導出するため、
    単純差なら solver の β_e も平均電流の中性も**厳密に 0** になり自己整合する
    （β≲0.2 では相対論補正 O(β²) は M_S の丸め誤差以下）。

    相対論的速度合成則が必要なほどドリフトが速い場合は relativistic=True。この場合
    フレーム速度（電子の速度）は厳密に 0 だが、solver が再導出する β_e には O(β²)
    の残差が残る（相対論では電流中性がフレーム依存になるため）。

    温度は非相対論ブーストでは不変とみなし、そのまま保つ。返り値に β_e キーは
    含めない（元々 solver 入力に β_e は無く、ソルバーが内部で再導出する）。

    Note
    ----
    平衡を電子静止系で**解き直す**には solve_snr_equilibrium_electron_frame を使う。
    これは自然系で解いた解を warm-start にするフレーム継続で β_e=0 を扱い、下記の
    seed 縮退を回避する。本関数の出力を **solve_snr_equilibrium に直接渡してはいけない**:
    線形分散 seed が結合係数 C_AP = β_e/Te + Σ_s w_s β_s/T_s に依存し、Tinc=Tref のとき
    C_AP = β_e·(1/Te + 1/T_i) なので β_e=0 で C_AP=0 となり、電磁モードの種が縮退して
    純静電の短波長分枝（誤った平衡）に落ちる。
    本関数は **線形解析（電子静止系での感受率）用にドリフトを準備する**もの。

    Parameters
    ----------
    params:
        build_solver_params の出力 {eta, beta_inc, beta_ref, Te, Tinc, Tref}。
    relativistic:
        False（既定）なら単純差 β - β_e。線形性により solver の β_e と平均電流の
        中性が厳密に保たれ、非相対論平衡と自己整合する。
        True なら相対論的速度合成則。電子のフレーム速度は厳密に 0 だが、平均電流の
        中性は O(β²) まで。

    Returns
    -------
    dict:
        electron-rest-frame に移したソルバー入力。eta/温度は不変、
        beta_inc/beta_ref のみブーストされる。
    """
    eta = params["eta"]
    beta_inc = params["beta_inc"]
    beta_ref = params["beta_ref"]
    beta_e = (1.0 - eta) * beta_inc + eta * beta_ref  # 電流中性条件

    out: dict[str, float] = dict(params)
    if relativistic:
        out["beta_inc"] = _boost_velocity(beta_inc, beta_e)
        out["beta_ref"] = _boost_velocity(beta_ref, beta_e)
    else:
        out["beta_inc"] = beta_inc - beta_e
        out["beta_ref"] = beta_ref - beta_e
    return out


def build_pic_config(cfg: MachConfig) -> dict:
    """MachConfig から PIC [parameter] 値を返す（新規ラン設計用）。

    磁化時(sigma>0)は betae/betai/betar（プラズマβ）で温度を表現できる。
    非磁化時(sigma=0)は v_A=0 で β が定義できないため、温度を直接指定する
    初期化（main.cpp 小改修）が必要 → vte/vti/vtr と needs_direct_temperature=True
    を返す。

    Parameters
    ----------
    cfg: マッハ数基準のパラメータ点。

    Returns
    -------
    dict:
        PIC config の [parameter] に書ける値。磁化時は betae/betai/betar、
        非磁化時は vte_over_c/vti_over_c/vtr_over_c を含む。
        後者のキー名は PIC 側変数（main.cpp の vte/vti/vtr）に対応する外部規約で、
        値は β_th=v_th/c（= 本モジュールの β_t* と同値）。
    """
    beta_sh = _beta_shock(cfg)
    beta_d = beta_sh / 2.0
    ush = beta_d / math.sqrt(1.0 - beta_d * beta_d)  # 4-velocity（β_d=ush/√(1+ush²) の逆）

    out: dict = dict(mime=cfg.mime, alpha=cfg.eta, sigma=cfg.sigma, ush=float(ush))
    sp = build_solver_params(cfg)
    # 熱速度（c 単位）: β_th,i=√(T̃_inc/mime), β_th,e=√T̃_e
    beta_ti = math.sqrt(sp["Tinc"] / cfg.mime)
    beta_tr = math.sqrt(sp["Tref"] / cfg.mime)
    beta_te = math.sqrt(sp["Te"])

    if cfg.sigma > 0.0:
        beta_ae2 = cfg.sigma  # (v_A,e/c)² = β_A,e²
        beta_ai2 = cfg.sigma / cfg.mime  # (v_A,i/c)² = β_A,i²
        out.update(
            betae=float(2.0 * beta_te ** 2 / beta_ae2),  # = 現行ラン 8.0 を再現
            betai=float(2.0 * beta_ti ** 2 / beta_ai2),  # = (2/γ)(M_A/M_S)² ≈ 0.5
            betar=float(2.0 * beta_tr ** 2 / beta_ai2),
            needs_direct_temperature=False,
        )
    else:
        # PIC 側変数名（main.cpp: vte/vti/vtr）の外部規約。値は β_th=v_th/c。
        out.update(
            vte_over_c=float(beta_te),
            vti_over_c=float(beta_ti),
            vtr_over_c=float(beta_tr),
            needs_direct_temperature=True,  # main.cpp の温度直接指定改修が必要
        )
    return out


__all__ = [
    "MachConfig",
    "build_solver_params",
    "build_pic_config",
    "to_electron_rest_frame",
]
