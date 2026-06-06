import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from fmi.snr_equilibrium import solve_snr_equilibrium
    return mo, np, plt, solve_snr_equilibrium


@app.cell
def _(mo):
    mo.md("""
    # SNR 3成分プラズマ 電流フィラメント平衡解

    Vanthieghem et al. (2018) のChebyshevスペクトル法 + Newton継続法を
    非相対論的3成分プラズマ（電子・入射イオン・反射イオン）に拡張。

    **正規化**: 長さ $c/\\omega_{pe}$、ポテンシャル $m_e c^2/e$、温度 $m_e c^2/k_B$
    """)
    return


@app.cell
def _(mo):
    eta_slider      = mo.ui.slider(0.05, 0.5,  step=0.05, value=0.1,  label="η（反射イオン割合）")
    beta_inc_slider = mo.ui.slider(0.05, 0.5,  step=0.05, value=0.3,  label="β_inc")
    beta_ref_slider = mo.ui.slider(-0.8, -0.1, step=0.05, value=-0.6, label="β_ref")
    Te_slider       = mo.ui.slider(0.1,  2.0,  step=0.1,  value=0.5,  label="T_e [m_e c²/k_B]")
    Tinc_slider     = mo.ui.slider(0.1,  2.0,  step=0.1,  value=1.0,  label="T_inc")
    Tref_slider     = mo.ui.slider(0.1,  2.0,  step=0.1,  value=0.5,  label="T_ref")
    a0_slider       = mo.ui.slider(0.05, 0.5,  step=0.05, value=0.25, label="a₀（目標振幅）")
    mo.vstack([
        mo.hstack([eta_slider, beta_inc_slider, beta_ref_slider]),
        mo.hstack([Te_slider, Tinc_slider, Tref_slider, a0_slider]),
    ])
    return eta_slider, beta_inc_slider, beta_ref_slider, Te_slider, Tinc_slider, Tref_slider, a0_slider


@app.cell
def _(mo, eta_slider, beta_inc_slider, beta_ref_slider,
      Te_slider, Tinc_slider, Tref_slider, a0_slider,
      solve_snr_equilibrium):
    with mo.status.spinner("平衡解を計算中..."):
        result = solve_snr_equilibrium(
            eta      = eta_slider.value,
            beta_inc = beta_inc_slider.value,
            beta_ref = beta_ref_slider.value,
            Te       = Te_slider.value,
            Tinc     = Tinc_slider.value,
            Tref     = Tref_slider.value,
            a0_target= a0_slider.value,
        )
    return result,


@app.cell
def _(mo, result):
    r = result
    mo.md(f"""
    **電流中性条件から導出**:
    $\\beta_e = (1-\\eta)\\beta_{{inc}} + \\eta\\beta_{{ref}} = {r['beta_e']:.4f}$
    &emsp;|&emsp;
    $\\lambda_0 = {r['lambda0']:.3f}\\; c/\\omega_{{pe}}$
    """)
    return


@app.cell
def _(plt, result):
    d = result
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # 左パネル: 電磁場
    ax1.plot(d["B0z_full"], d["y_full"], "k-",  lw=1.5, label=r"$B_{0z}$")
    ax1.plot(d["E0y_full"], d["y_full"], "k--", lw=1.5, label=r"$E_{0y}$")
    ax1.set_xlabel(r"Fields $[m_e c\,\omega_{pe}/e]$")
    ax1.set_ylabel(r"$y\;[c/\omega_{pe}]$")
    ax1.set_title("Electromagnetic Fields")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 右パネル: 密度プロファイル（n̄_e で正規化）
    ax2.plot(d["ne_full"],   d["y_full"], "b-",  lw=1.5, label=r"$n_e / \bar{n}_e$")
    ax2.plot(d["ninc_full"], d["y_full"], "r-",  lw=1.5, label=r"$n_{inc} / \bar{n}_e$")
    ax2.plot(d["nref_full"], d["y_full"], "g-",  lw=1.5, label=r"$n_{ref} / \bar{n}_e$")
    ax2.set_xlabel(r"Density $[\bar{n}_e]$")
    ax2.set_title("Density Profiles")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    p = d["params"]
    fname = (f"work/snr_eq_eta{p['eta']:.2f}_binc{p['beta_inc']:.2f}"
             f"_bref{p['beta_ref']:.2f}_a0{p['a0_target']:.2f}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    fig
    return fig,


if __name__ == "__main__":
    app.run()
