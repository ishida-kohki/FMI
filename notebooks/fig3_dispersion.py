import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from fmi import equilibrium, floquet

    return equilibrium, floquet, mo, np, plt


@app.cell
def _(mo):
    mo.md(
        """
        # FMI 分散関係 — Fig. 3 (Vanthieghem et al. 2018)

        Floquet 解析により、純粋成長モード (omega = i Gamma Omega_p, kx=0) の
        分散関係 Gamma(ky) を計算します。
        """
    )
    return


@app.cell
def _(mo):
    xi_slider = mo.ui.slider(0.4, 5.0, step=0.1, value=5.0, label="xi")
    n_scan_slider = mo.ui.slider(20, 120, step=10, value=80, label="Gamma scan points")
    mo.hstack([xi_slider, n_scan_slider])
    return n_scan_slider, xi_slider


@app.cell
def _(equilibrium, n_scan_slider, np, xi_slider):
    xi = xi_slider.value

    params, Omega_p, k0_paper = equilibrium.build_equilibrium(xi=xi)

    gamma_max_guess = (2.0 / np.sqrt(5.0)) * np.exp(-xi / 2.0)
    gamma_paper_list = np.linspace(gamma_max_guess * 1.1, 0.001, n_scan_slider.value)

    return Omega_p, gamma_paper_list, k0_paper, params, xi


@app.cell
def _(Omega_p, floquet, gamma_paper_list, k0_paper, mo, params):
    with mo.status.spinner("Floquet scan 実行中..."):
        ky_plot, gamma_plot = floquet.scan_dispersion(
            kx=0.0,
            params=params,
            Omega_p=Omega_p,
            k0_paper=k0_paper,
            gamma_list=gamma_paper_list,
        )

    return gamma_plot, ky_plot


@app.cell
def _(gamma_plot, k0_paper, ky_plot, plt, xi):
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(ky_plot, gamma_plot, s=12, color='blue', label=rf'FMI ($\xi={xi}$)')
    ax.axvline(k0_paper / 2.0, color='red',   linestyle='--', linewidth=1.2, alpha=0.7, label=r'$k_0/2$ (Gap)')
    ax.axvline(k0_paper,       color='green',  linestyle='--', linewidth=1.2, alpha=0.7, label=r'$k_0$')

    ax.set_xlabel(r'Transverse wavenumber $k_y\ [c/\Omega_p]$', fontsize=13)
    ax.set_ylabel(r'Growth rate $\Gamma\ [\Omega_p]$', fontsize=13)
    ax.set_title(rf'Dispersion Relation (Floquet Analysis, $\xi={xi}$)', fontsize=14)
    ax.set_xlim(0, k0_paper * 1.1)
    ax.set_ylim(0, None)
    ax.grid(True, linestyle=':')
    ax.legend(loc='upper right', fontsize=12)
    fig.tight_layout()
    fig.savefig(f"work/fig3_xi{xi}.png", dpi=300, bbox_inches='tight')
    fig
    return (fig,)


@app.cell
def _(Omega_p, k0_paper, mo, params, xi):
    mo.md(
        f"""
        **パラメータ**

        | 変数 | 値 |
        |------|-----|
        | xi | {xi} |
        | lambda0 | {params["lambda0"]:.4f} |
        | Omega_p | {Omega_p:.4f} |
        | k0 (normalized) | {k0_paper:.4f} Omega_p/c |
        | k0/2 (gap) | {k0_paper/2:.4f} Omega_p/c |
        """
    )
    return


if __name__ == "__main__":
    app.run()
