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
        # FMI 分散関係 — Fig. 10 (Vanthieghem et al. 2018)

        複数の kx 値における分散関係 Gamma(ky) を Floquet 解析で計算します。
        """
    )
    return


@app.cell
def _(mo):
    xi_slider = mo.ui.slider(0.4, 5.0, step=0.1, value=5.0, label="xi")
    n_scan_slider = mo.ui.slider(20, 120, step=10, value=60, label="Gamma scan points")
    mo.hstack([xi_slider, n_scan_slider])
    return n_scan_slider, xi_slider


@app.cell
def _(equilibrium, n_scan_slider, np, xi_slider):
    xi = xi_slider.value

    params, Omega_p, k0_paper = equilibrium.build_equilibrium(xi=xi)

    gamma_list = np.linspace(0.002, 0.25, n_scan_slider.value)
    kx_norm_list = [0.0, 0.03, 0.16, 0.42]

    return Omega_p, gamma_list, k0_paper, kx_norm_list, params, xi


@app.cell
def _(Omega_p, floquet, gamma_list, k0_paper, kx_norm_list, mo, params):
    results = {}

    with mo.status.spinner("Floquet scan 実行中 (複数 kx)..."):
        for kx_norm in kx_norm_list:
            ky_plot, gamma_plot = floquet.scan_dispersion(
                kx=kx_norm * Omega_p,
                params=params,
                Omega_p=Omega_p,
                k0_paper=k0_paper,
                gamma_list=gamma_list,
            )
            results[kx_norm] = (ky_plot, gamma_plot)

    return (results,)


@app.cell
def _(k0_paper, kx_norm_list, plt, results, xi):
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ['C0', 'C1', 'C2', 'C3']

    for i, kx_norm in enumerate(kx_norm_list):
        ky_plot, gamma_plot = results[kx_norm]
        if ky_plot:
            ax.scatter(ky_plot, gamma_plot, s=5, c=colors[i % len(colors)],
                       alpha=0.5, label=rf'$k_x={kx_norm}\,\Omega_p/c$')

    ax.axvline(k0_paper / 2.0, color='k',    linestyle='--', lw=1.2, label=r'$k_0/2$')
    ax.axvline(k0_paper,       color='gray',  linestyle=':',  lw=1.2, label=r'$k_0$')

    ax.set_xlabel(r'$k_y\ [\Omega_p/c]$', fontsize=13)
    ax.set_ylabel(r'$\Gamma\ [\Omega_p]$', fontsize=13)
    ax.set_title(rf'Dispersion relation — Floquet, $\xi={xi}$, multiple $k_x$', fontsize=13)
    ax.set_xlim(0.0, k0_paper * 1.1)
    ax.set_ylim(0.0, 0.25)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig
    return (fig,)


if __name__ == "__main__":
    app.run()
