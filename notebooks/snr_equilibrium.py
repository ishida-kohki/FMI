import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from fmi.snr_equilibrium import (
        solve_snr_equilibrium,
        solve_snr_equilibrium_electron_frame,
    )

    return mo, plt, solve_snr_equilibrium, solve_snr_equilibrium_electron_frame


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # SNR 3-Component Plasma Current Filament Equilibrium

    Chebyshev spectral method + Newton continuation (Vanthieghem et al. 2018)
    extended to a non-relativistic 3-component plasma (electrons, incident ions, reflected ions).

    **Units**: length $c/\omega_{pe}$, potential $m_e c^2/e$, temperature $m_e c^2/k_B$
    """)
    return


@app.cell(hide_code=True)
def _(mo, pp):
    pp
    # 温度・ドリフトは number 入力（任意符号・任意精度）にして
    # マッハ数由来の値（T~1e-4, β_inc<0, β_ref>0 など）を厳密に入力できるようにする。
    eta_slider      = mo.ui.slider(0.05, 0.9, step=0.05, value=0.2, label="η (reflected ion fraction)")
    beta_inc_slider = mo.ui.number(-1.0, 1.0, value=0.3,  label="β_inc")
    beta_ref_slider = mo.ui.number(-1.0, 1.0, value=-0.6, label="β_ref")
    # T は 1e-5 まで入力可（step 省略 = 任意値）。マッハ既定の T_inc≈6.25e-4 等を直接入力。
    Te_slider       = mo.ui.number(1e-5, 2.0, value=0.5,  label="T_e [m_e c²/k_B]")
    Tinc_slider     = mo.ui.number(1e-5, 2.0, value=1.0,  label="T_inc")
    Tref_slider     = mo.ui.number(1e-5, 2.0, value=0.5,  label="T_ref")
    # a₀ は数値入力（0.01–100）。高 a₀ ほど非線形 ξ=βa₀/T が上がるが、ξ≳12 でシートが
    # 格子より薄くなり under-resolved、ξ>100 で exp クリップ発動で解が破綻する点に注意。
    a0_slider       = mo.ui.number(0.01, 100.0, value=0.25, label="a₀ (target amplitude)")
    frame_dd        = mo.ui.dropdown(
        ["natural (β_e≠0)", "electron rest (β_e=0)"],
        value="natural (β_e≠0)", label="solver frame",
    )
    mo.vstack([
        mo.hstack([eta_slider, beta_inc_slider, beta_ref_slider]),
        mo.hstack([Te_slider, Tinc_slider, Tref_slider, a0_slider]),
        mo.hstack([frame_dd]),
    ])
    return (
        Te_slider,
        Tinc_slider,
        Tref_slider,
        a0_slider,
        beta_inc_slider,
        beta_ref_slider,
        eta_slider,
        frame_dd,
    )


@app.cell(hide_code=True)
def _(mo):
    # マッハ既定（M_S=154.92, η=0.2, V_sh/c=0.25, v_te/c=0.1, σ=0）の自然系入力値を表示。
    # 上の入力欄にこれをそのまま入れ、frame=electron rest を選べばマッハ図を再現できる。
    from fmi.mach_parameters import MachConfig, build_solver_params
    _sp = build_solver_params(
        MachConfig(M_S=154.92, beta_sh=0.25, beta_te=0.1,
                   sigma=0.0, eta=0.2, mime=400.0, gamma=5.0 / 3.0)
    )
    mo.md(f"""
    **マッハ由来の入力値（自然系, η=0.2）** — 上の欄にこれを入力して比較:
    β_inc = `{_sp['beta_inc']:.6f}`, β_ref = `{_sp['beta_ref']:.6f}`,
    T_e = `{_sp['Te']:.6e}`, T_inc = `{_sp['Tinc']:.6e}`, T_ref = `{_sp['Tref']:.6e}`
    &emsp;（マッハ図は **electron rest** フレームで生成。frame を合わせること）
    """)
    return


@app.cell(hide_code=True)
def _(
    Te_slider,
    Tinc_slider,
    Tref_slider,
    a0_slider,
    beta_inc_slider,
    beta_ref_slider,
    eta_slider,
    frame_dd,
    mo,
    solve_snr_equilibrium,
    solve_snr_equilibrium_electron_frame,
):
    _solver = (solve_snr_equilibrium_electron_frame
               if frame_dd.value.startswith("electron")
               else solve_snr_equilibrium)
    with mo.status.spinner(f"Computing equilibrium ({frame_dd.value})..."):
        result = _solver(
            eta      = eta_slider.value,
            beta_inc = beta_inc_slider.value,
            beta_ref = beta_ref_slider.value,
            Te       = Te_slider.value,
            Tinc     = Tinc_slider.value,
            Tref     = Tref_slider.value,
            a0_target= a0_slider.value,
        )
    return (result,)


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo, result):
    # 変数は _ 接頭辞でセルローカルにする（marimo は同名のグローバル再定義を禁止）
    _r = result
    _p = _r["params"]
    _xi_inc = abs(_p["beta_inc"]) * _p["a0_target"] / _p["Tinc"]
    _xi_ref = abs(_p["beta_ref"]) * _p["a0_target"] / _p["Tref"]
    _eb = abs(_r["E0y_full"]).max() / max(abs(_r["B0z_full"]).max(), 1e-30)
    mo.md(f"""
    **Derived from current neutrality**:
    $\\beta_e = (1-\\eta)\\beta_{{inc}} + \\eta\\beta_{{ref}} = {_r['beta_e']:.4f}$
    &emsp;|&emsp;
    $\\lambda_0 = {_r['lambda0']:.3f}\\; c/\\omega_{{pe}}$

    **Nonlinearity** $\\xi_s = |\\beta_s|\\,a_0/T_s$ (Harris は $\\xi\\gg1$):
    $\\xi_{{inc}} = {_xi_inc:.2f}$ &emsp; $\\xi_{{ref}} = {_xi_ref:.2f}$
    &emsp;|&emsp;
    $|E_{{0y}}|_{{max}}/|B_{{0z}}|_{{max}} = {_eb:.3f}$ (対称なら ambipolar 場のみ)
    """)
    return


@app.cell(hide_code=True)
def _(plt, result):
    d = result
    p = d["params"]

    # 電荷密度 ρ と電流密度 j_x（ソルバ正規化と一致, q_e=-1, q_ion=+1）
    #   ρ   = -n_e + n_inc + n_ref
    #   j_x = -β_e n_e + β_inc n_inc + β_ref n_ref   （dB_0z/dy = j_x）
    rho = -d["ne_full"] + d["ninc_full"] + d["nref_full"]
    jx = (-d["beta_e"] * d["ne_full"]
          + p["beta_inc"] * d["ninc_full"]
          + p["beta_ref"] * d["nref_full"])
    n_ion = d["ninc_full"] + d["nref_full"]   # 総イオン密度（Harris sech² 比較の量）

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 14), sharex=True)

    ax1.plot(d["y_full"], d["B0z_full"], "k-",  lw=1.5, label=r"$B_{0z}$")
    ax1.plot(d["y_full"], d["E0y_full"], "k--", lw=1.5, label=r"$E_{0y}$")
    ax1.axhline(0, color="gray", lw=0.5, ls=":")
    ax1.set_ylabel(r"Fields $[m_e c\,\omega_{pe}/e]$")
    ax1.set_title("Electromagnetic Fields")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(d["y_full"], d["ne_full"],   "b-",  lw=1.5, label=r"$n_e / \bar{n}_e$")
    ax2.plot(d["y_full"], d["ninc_full"], "r-",  lw=1.5, label=r"$n_{inc} / \bar{n}_e$")
    ax2.plot(d["y_full"], d["nref_full"], "g-",  lw=1.5, label=r"$n_{ref} / \bar{n}_e$")
    ax2.set_ylabel(r"Density $[\bar{n}_e]$")
    ax2.set_title("Density Profiles")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    ax3.plot(d["y_full"], n_ion, "m-", lw=1.5, label=r"$n_{ion} = n_{inc} + n_{ref}$")
    ax3.set_ylabel(r"$n_{ion}\;[\bar{n}_e]$")
    ax3.set_title(r"Total Ion Density (Harris $\mathrm{sech}^2$ comparison)")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    ax4.plot(d["y_full"], rho, "k-",  lw=1.5, label=r"$\rho = -n_e + n_{inc} + n_{ref}$")
    ax4.plot(d["y_full"], jx,  "k--", lw=1.5, label=r"$j_x = -\beta_e n_e + \beta_{inc} n_{inc} + \beta_{ref} n_{ref}$")
    ax4.axhline(0, color="gray", lw=0.5, ls=":")
    ax4.set_xlabel(r"$y\;[c/\omega_{pe}]$")
    ax4.set_ylabel(r"$\rho\;[e\bar{n}_e]\;,\;\;j_x\;[e\bar{n}_e c]$")
    ax4.set_title("Charge & Current Density")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    fig.tight_layout()
    fname = (f"work/snr_eq_eta{p['eta']:.2f}_binc{p['beta_inc']:.2f}"
             f"_bref{p['beta_ref']:.2f}_a0{p['a0_target']:.2f}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    fig
    return


if __name__ == "__main__":
    app.run()
