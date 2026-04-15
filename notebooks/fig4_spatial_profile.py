import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import scipy.linalg as la
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    from fmi import equilibrium
    from fmi.floquet import get_N_matrix

    return equilibrium, get_N_matrix, la, mo, np, plt, solve_ivp


@app.cell
def _(mo):
    mo.md(
        """
        # FMI 空間プロファイル — Fig. 4 (Vanthieghem et al. 2018)

        最大成長率をもつ FMI モードの空間プロファイル
        (δEx, δEy, δBz, δv₁y, δd₁) を Floquet 固有ベクトルから再構築します。
        """
    )
    return


@app.cell
def _(mo):
    xi_slider     = mo.ui.slider(0.4, 5.0, step=0.1,  value=0.4,  label="xi")
    gamma_slider  = mo.ui.slider(0.1, 1.5, step=0.01, value=0.685, label="Gamma / Omega_p")
    y_max_slider  = mo.ui.slider(1.0, 8.0, step=0.1,  value=3.8,  label="y_max [c/Omega_p]")
    mo.hstack([xi_slider, gamma_slider, y_max_slider])
    return gamma_slider, xi_slider, y_max_slider


@app.cell
def _(equilibrium, gamma_slider, get_N_matrix, la, mo, np, solve_ivp, xi_slider, y_max_slider):
    xi    = xi_slider.value
    params, Omega_p, k0_paper = equilibrium.build_equilibrium(xi=xi)
    lambda0_native = params['lambda0']
    kx_test = 0.0

    Gamma_target = gamma_slider.value * Omega_p
    omega_test   = 1j * Gamma_target

    # --- 固有ベクトル抽出 (Floquet モノドロミー行列) ---
    X_init_r = np.concatenate([np.eye(10, dtype=np.complex128).real.flatten(),
                                np.eye(10, dtype=np.complex128).imag.flatten()])

    def ode_system_real(y, X_flat_real):
        X_R = X_flat_real[:100].reshape((10, 10))
        X_I = X_flat_real[100:].reshape((10, 10))
        N = get_N_matrix(y, omega_test, kx_test, params)
        dX_R = N.real @ X_R - N.imag @ X_I
        dX_I = N.real @ X_I + N.imag @ X_R
        return np.concatenate([dX_R.flatten(), dX_I.flatten()])

    with mo.status.spinner("Floquet 固有ベクトルを計算中..."):
        sol_fwd = solve_ivp(ode_system_real, [lambda0_native / 2.0, lambda0_native],
                            X_init_r, method='Radau', rtol=1e-8, atol=1e-8)
        sol_bwd = solve_ivp(ode_system_real, [lambda0_native / 2.0, 0.0],
                            X_init_r, method='Radau', rtol=1e-8, atol=1e-8)

    X_L = sol_fwd.y[:100, -1].reshape((10, 10)) + 1j * sol_fwd.y[100:, -1].reshape((10, 10))
    X_0 = sol_bwd.y[:100, -1].reshape((10, 10)) + 1j * sol_bwd.y[100:, -1].reshape((10, 10))
    B = la.inv(X_0) @ X_L
    eigvals, eigvecs = la.eig(B)

    dist_from_unit = np.abs(np.abs(eigvals) - 1.0)
    valid_idx = np.where(dist_from_unit < 0.15)[0]
    if len(valid_idx) == 0:
        best_idx = np.argmin(dist_from_unit)
    else:
        angles = np.angle(eigvals[valid_idx])
        pos_mask = angles > 0.1
        pos_idx = valid_idx[pos_mask]
        if len(pos_idx) > 0:
            best_idx = pos_idx[np.argmax(np.angle(eigvals[pos_idx]))]
        else:
            best_idx = valid_idx[np.argmax(angles)]

    best_eigvec = eigvecs[:, best_idx]
    V_y0 = X_0 @ best_eigvec
    # Re(δBz) がピーク(+1)になるように位相を揃える
    V_y0_aligned = V_y0 * np.exp(-1j * np.angle(V_y0[1]))

    # --- 1周期の空間積分 (dense_output=True でタイリング用) ---
    def ode_single_mode(y, X_flat):
        X_c = X_flat[:10] + 1j * X_flat[10:]
        N = get_N_matrix(y, omega_test, kx_test, params)
        dX_c = N @ X_c
        return np.concatenate([dX_c.real, dX_c.imag])

    V_r = np.concatenate([V_y0_aligned.real, V_y0_aligned.imag])

    with mo.status.spinner("空間プロファイルを積分中..."):
        sol_piece = solve_ivp(ode_single_mode, [0.0, lambda0_native], V_r,
                              dense_output=True, method='Radau', rtol=1e-8, atol=1e-8)

    mu_raw    = eigvals[best_idx]
    mu_phase  = mu_raw / np.abs(mu_raw)
    decay_rate = -np.log(np.abs(mu_raw)) / lambda0_native

    y_max_native = y_max_slider.value / Omega_p
    t_eval = np.linspace(0, y_max_native, 1000)
    y_plot = t_eval * Omega_p  # [c/Omega_p] 単位

    X_sol_c = np.zeros((10, len(t_eval)), dtype=np.complex128)
    for i, y_val in enumerate(t_eval):
        n_periods = int(np.floor(y_val / lambda0_native))
        y_mod = np.clip(y_val % lambda0_native, sol_piece.t[0], sol_piece.t[-1])
        X_flat_mod = sol_piece.sol(y_mod)
        X_c_mod = X_flat_mod[:10] + 1j * X_flat_mod[10:]
        X_c_mod_corrected = X_c_mod * np.exp(decay_rate * y_mod)
        X_sol_c[:, i] = X_c_mod_corrected * (mu_phase ** n_periods)

    # --- 場の抽出 (index: 0=δEx, 1=δBz, 2+a=δPa, 6+a=δVay) ---
    delta_Ex  = X_sol_c[0, :]
    delta_Bz  = X_sol_c[1, :]
    delta_P1  = X_sol_c[3, :]   # a=1
    delta_V1y = X_sol_c[7, :]   # a=1

    # δEy と δd₁ を代数関係式から再構築
    bg_func  = params['bg_func']
    q_arr    = params['q']
    eps_arr  = params['eps']
    beta0    = params['beta0']
    gamma0   = params['gamma0']
    d0       = params['d0']
    p0       = params['p0']
    G_ad     = params['Gamma_ad']
    c        = params['c']
    m        = params['m']

    delta_Ey = np.zeros_like(delta_Ex, dtype=np.complex128)
    delta_d1 = np.zeros_like(delta_Ex, dtype=np.complex128)
    bg_B0z   = np.zeros(len(t_eval))
    bg_d1    = np.zeros(len(t_eval))

    for i, y_val in enumerate(t_eval):
        B0z_c, _, f_B, _, f_a = bg_func(y_val)
        bg_B0z[i] = B0z_c * f_B
        bg_d1[i]  = f_a[1]

        sum_Vy = sum(q_arr[b] * beta0[b] * d0[b] * f_a[b] * X_sol_c[6 + b, i] for b in range(4))
        delta_Ey[i] = -1j / omega_test * (4.0 * np.pi * c / B0z_c) * sum_Vy

        h0_1 = p0[1] * G_ad[1] / (G_ad[1] - 1.0) + m * c**2 * d0[1] / gamma0[1]
        W1   = (gamma0[1]**2 / c**2) * h0_1 * f_a[1] * beta0[1] * (-1j * omega_test)
        C_Vx_ex_1  = q_arr[1] * d0[1] * B0z_c * f_a[1] * (1.0 - beta0[1]**2) / W1
        C_Vx_P_1   = (1j * eps_arr[1] * omega_test * p0[1] * beta0[1] / c) / W1
        C_Vx_Vy_1  = q_arr[1] * d0[1] * beta0[1] * f_a[1] * (B0z_c * f_B) / W1
        K_1        = eps_arr[1] * f_a[1] * beta0[1]**2 / (1.0 - beta0[1]**2)
        delta_d1[i] = (1.0 / G_ad[1]) * delta_P1[i] + K_1 * (
            C_Vx_ex_1 * delta_Ex[i] + C_Vx_P_1 * delta_P1[i] + C_Vx_Vy_1 * delta_V1y[i]
        )

    # --- 規格化 ---
    max_Bz = np.max(np.abs(delta_Bz))
    delta_Ex /= max_Bz;  delta_Ey /= max_Bz;  delta_Bz /= max_Bz
    max_d1 = np.max(np.abs(delta_d1))
    delta_V1y /= max_d1; delta_d1 /= max_d1
    bg_B0z /= np.max(np.abs(bg_B0z))
    bg_d1  /= np.max(np.abs(bg_d1))

    return (
        bg_B0z, bg_d1, delta_Bz, delta_Ex, delta_Ey, delta_V1y, delta_d1,
        xi, y_plot,
    )


@app.cell
def _(bg_B0z, bg_d1, delta_Bz, delta_Ex, delta_Ey, delta_V1y, delta_d1, gamma_slider, plt, xi, y_max_slider, y_plot):
    fig, axs = plt.subplots(1, 4, figsize=(10, 6), sharey=True)

    axs[0].plot(delta_Ex.real,  y_plot, 'k-', lw=2,        label=r'Re($\delta E_x$)')
    axs[0].plot(delta_Ex.imag,  y_plot, 'k-', lw=1, alpha=0.5, label=r'Im($\delta E_x$)')
    axs[0].plot(delta_Ey.real,  y_plot, 'r-', lw=2,        label=r'Re($\delta E_y$)')
    axs[0].plot(delta_Ey.imag,  y_plot, 'r-', lw=1, alpha=0.5, label=r'Im($\delta E_y$)')
    axs[0].axvline(0, color='gray', lw=1)
    axs[0].set_xlabel(r'$\delta E$')
    axs[0].set_ylabel(r'$y\ [c/\Omega_p]$')
    axs[0].legend(loc='lower right', fontsize=7)

    axs[1].plot(delta_Bz.real,  y_plot, 'k-',  lw=2,        label=r'Re($\delta B_z$)')
    axs[1].plot(delta_Bz.imag,  y_plot, 'k-',  lw=1, alpha=0.5, label=r'Im($\delta B_z$)')
    axs[1].plot(bg_B0z,         y_plot, 'b--', lw=1, alpha=0.6, label=r'$B_{0z}$ (unperturbed)')
    axs[1].axvline(0, color='gray', lw=1)
    axs[1].set_xlabel(r'$\delta B_z$')
    axs[1].legend(loc='lower right', fontsize=7)

    axs[2].plot(delta_V1y.real, y_plot, 'k-', lw=2,        label=r'Re($\delta v_{1y}$)')
    axs[2].plot(delta_V1y.imag, y_plot, 'k-', lw=1, alpha=0.5, label=r'Im($\delta v_{1y}$)')
    axs[2].axvline(0, color='gray', lw=1)
    axs[2].set_xlabel(r'$\delta v_{1y}$')
    axs[2].legend(loc='lower right', fontsize=7)

    axs[3].plot(delta_d1.real,  y_plot, 'k-',  lw=2,        label=r'Re($\delta d_1$)')
    axs[3].plot(delta_d1.imag,  y_plot, 'k-',  lw=1, alpha=0.5, label=r'Im($\delta d_1$)')
    axs[3].plot(bg_d1,          y_plot, 'b--', lw=1, alpha=0.6, label=r'$d_{1,0}$ (unperturbed)')
    axs[3].axvline(0, color='gray', lw=1)
    axs[3].set_xlabel(r'$\delta d_1$')
    axs[3].legend(loc='lower right', fontsize=7)

    for ax in axs:
        ax.set_ylim(0, y_max_slider.value)
        ax.tick_params(direction='in')

    fig.suptitle(rf'Spatial profile — $\xi={xi}$, $\Gamma={gamma_slider.value}\,\Omega_p$', fontsize=12)
    fig.tight_layout()
    fig.savefig(f"work/fig4_xi{xi}_gamma{gamma_slider.value:.3f}.png", dpi=300, bbox_inches='tight')
    fig
    return (fig,)


if __name__ == "__main__":
    app.run()
