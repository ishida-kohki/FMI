import marimo

__generated_with = "0.20.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    from scipy.optimize import root
    import matplotlib.pyplot as plt

    return mo, np, plt, root


@app.cell
def _(mo):
    mo.md(
        """
        # FMI 均衡プロファイル — Fig. 1 (Vanthieghem et al. 2018)

        Chebyshev スペクトル法 + Newton 継続法により、非対称均衡
        （ビームプラズマ系）の電磁場および密度プロファイルを求めます。
        """
    )
    return


@app.cell
def _(mo):
    a0_slider = mo.ui.slider(0.05, 0.5, step=0.05, value=0.25, label="a0 target")
    mo.hstack([a0_slider])
    return (a0_slider,)


@app.cell
def _(a0_slider, mo, np, root):
    def chebyshev_diff_matrix(N):
        if N == 0:
            return np.array([[0.0]]), np.array([0.0])
        x = np.cos(np.arange(N + 1) * np.pi / N)
        D = np.zeros((N + 1, N + 1))
        p = np.ones(N + 1); p[0] = 2; p[N] = 2
        for i in range(N + 1):
            for j in range(N + 1):
                if i != j:
                    D[i, j] = (-1)**(i + j) * p[i] / (p[j] * (x[i] - x[j]))
        D[0, 0] = (1 + 2 * N**2) / 6
        D[N, N] = -(1 + 2 * N**2) / 6
        for j in range(1, N):
            D[j, j] = -x[j] / (2 * (1 - x[j]**2))
        return D, x

    def solve_asymmetric_equilibrium(a0_target=0.25):
        N_points = 64
        """
        fig1の設定
        Tb = 1.0;  betab = -0.995; gammab = 1.0 / np.sqrt(1 - betab**2)
        Tp = 0.1;  betap =  0.995; gammap = 1.0 / np.sqrt(1 - betap**2)
        """
        Tb = 1.0;  betab = -0.995; gammab = 1.0 / np.sqrt(1 - betab**2)
        Tp = 0.1;  betap =  0.8; gammap = 1.0 / np.sqrt(1 - betap**2)
        D_std, x_std = chebyshev_diff_matrix(N_points)
        D2_std = D_std @ D_std

        def residual(u, current_a0):
            A   = u[0 : N_points+1]
            Phi = u[N_points+1 : 2*N_points+2]
            L_half = np.abs(u[-1]) + 1e-6
            dy_dx   = -(L_half / 2.0)
            D_phys  = D_std  / dy_dx
            D2_phys = D2_std / (dy_dx**2)
            arg_b = np.clip((gammab / Tb) * (Phi - betab * A), -100, 100)
            arg_p = np.clip((gammap / Tp) * (Phi - betap * A), -100, 100)
            rhs_A   = 2 * (gammab * betab * np.sinh(arg_b) + gammap * betap * np.sinh(arg_p))
            rhs_Phi = 2 * (gammab * np.sinh(arg_b)         + gammap * np.sinh(arg_p))
            eq_A   = (D2_phys @ A)   - rhs_A
            eq_Phi = (D2_phys @ Phi) - rhs_Phi
            factor = dy_dx**2
            res = np.zeros_like(u)
            res[1:N_points]                = eq_A[1:N_points]   * factor
            res[N_points+2 : 2*N_points+1] = eq_Phi[1:N_points] * factor
            res[0]            = (A[0] - current_a0)       * factor
            res[N_points+1]   = (D_phys @ Phi)[0]         * factor
            res[-1]           = (D_phys @ A)[0]           * factor
            res[N_points]     = (D_phys @ A)[N_points]    * factor
            res[2*N_points+1] = (D_phys @ Phi)[N_points]  * factor
            return res

        a0_start = 0.01
        C_AA     = -2 * (gammab**2 * betab**2 / Tb + gammap**2 * betap**2 / Tp)
        C_APhi   =  2 * (gammab**2 * betab   / Tb + gammap**2 * betap   / Tp)
        C_PhiPhi =  2 * (gammab**2           / Tb + gammap**2           / Tp)
        a_q = C_APhi; b_q = C_AA - C_PhiPhi; c_q = C_APhi
        D_q = b_q**2 - 4 * a_q * c_q
        c1 = (-b_q + np.sqrt(D_q)) / (2 * a_q)
        c2 = (-b_q - np.sqrt(D_q)) / (2 * a_q)
        if -(C_AA + c1 * C_APhi) > 0:
            phi_ratio = c1; k2 = -(C_AA + c1 * C_APhi)
        else:
            phi_ratio = c2; k2 = -(C_AA + c2 * C_APhi)
        L_half_guess = np.pi / np.sqrt(k2)
        u_guess = np.zeros(2 * N_points + 3)
        u_guess[0 : N_points+1]           = a0_start * np.cos(np.pi * (1 - x_std) / 2)
        u_guess[N_points+1 : 2*N_points+2] = phi_ratio * u_guess[0 : N_points+1]
        u_guess[-1] = L_half_guess

        for i, current_a0 in enumerate(np.linspace(a0_start, a0_target, 25)):
            if i > 0:
                u_guess[0 : 2*N_points+2] *= current_a0 / np.linspace(a0_start, a0_target, 25)[i-1]
            sol = root(residual, u_guess, args=(current_a0,), method='lm')
            u_guess = sol.x

        A   = sol.x[0 : N_points+1]
        Phi = sol.x[N_points+1 : 2*N_points+2]
        L_half = np.abs(sol.x[-1])
        dy_dx  = -(L_half / 2.0)
        D_phys = D_std / dy_dx
        B0z_hat = -(D_phys @ A)
        E0y_hat = -(D_phys @ Phi)
        arg_b = (gammab / Tb) * (Phi - betab * A)
        arg_p = (gammap / Tp) * (Phi - betap * A)
        np_minus = np.exp(arg_p); np_plus = np.exp(-arg_p)
        max_np = max(np.max(np_minus), np.max(np_plus))
        lambda0 = 2 * L_half * np.sqrt(max_np)
        B0z = B0z_hat / np.sqrt(max_np)
        E0y = E0y_hat / np.sqrt(max_np)
        y_half = L_half * (1 - x_std) / 2 * np.sqrt(max_np)
        db_minus = gammab * np.exp(arg_b)  / (gammab * np.exp(np.max(np.abs(arg_p))))
        db_plus  = gammab * np.exp(-arg_b) / (gammab * np.exp(np.max(np.abs(arg_p))))
        dp_minus = gammap * np.exp(arg_p)
        dp_plus  = gammap * np.exp(-arg_p)
        max_dp = max(np.max(dp_minus), np.max(dp_plus))
        db_minus /= max_dp / gammap; db_plus /= max_dp / gammap
        dp_minus_norm = dp_minus / max_dp; dp_plus_norm = dp_plus / max_dp
        db_minus_norm = gammab * np.exp(arg_b)  / max_dp
        db_plus_norm  = gammab * np.exp(-arg_b) / max_dp
        y_full        = np.concatenate([y_half, 2 * y_half[-1] - y_half[-2::-1]])
        B0z_full      = np.concatenate([B0z,  -B0z[-2::-1]])
        E0y_full      = np.concatenate([E0y,  -E0y[-2::-1]])
        db_minus_full = np.concatenate([db_minus_norm, db_minus_norm[-2::-1]])
        db_plus_full  = np.concatenate([db_plus_norm,  db_plus_norm[-2::-1]])
        dp_minus_full = np.concatenate([dp_minus_norm, dp_minus_norm[-2::-1]])
        dp_plus_full  = np.concatenate([dp_plus_norm,  dp_plus_norm[-2::-1]])
        return dict(
            y_full=y_full, B0z_full=B0z_full, E0y_full=E0y_full,
            db_minus_full=db_minus_full, db_plus_full=db_plus_full,
            dp_minus_full=dp_minus_full, dp_plus_full=dp_plus_full,
            lambda0=lambda0,
            Tb=Tb, Tp=Tp, betab=betab, betap=betap,
        )

    with mo.status.spinner("均衡解を計算中..."):
        result = solve_asymmetric_equilibrium(a0_target=a0_slider.value)

    return chebyshev_diff_matrix, result, solve_asymmetric_equilibrium


@app.cell
def _(a0_slider, plt, result):
    d = result
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharey=True)

    ax1.plot(d['B0z_full'], d['y_full'], 'k-',  label=r'$B_{0z}$')
    ax1.plot(d['E0y_full'], d['y_full'], 'k--', label=r'$E_{0y}$')
    ax1.set_title('Electromagnetic Fields')
    ax1.set_xlabel(r'Fields $\left[mc\Omega_p/e\right]$')
    ax1.set_ylabel(r'$y \ [c/\Omega_p]$')
    ax1.grid(True)
    ax1.legend(loc='lower right')

    ax2.plot(d['dp_minus_full'], d['y_full'], 'r-',  linewidth=1.5, label=r'$d_{p}^{-}$')
    ax2.plot(d['dp_plus_full'],  d['y_full'], 'r--', linewidth=1.5, label=r'$d_{p}^{+}$')
    ax2.plot(d['db_minus_full'], d['y_full'], 'b-',  linewidth=1.5, label=r'$d_{b}^{-}$')
    ax2.plot(d['db_plus_full'],  d['y_full'], 'b--', linewidth=1.5, label=r'$d_{b}^{+}$')
    ax2.set_title('Apparent Densities')
    ax2.set_xlabel(r'Density $d_{\alpha,0} / \max(d_{p,0})$')
    ax2.grid(True)
    ax2.legend(loc='lower right')

    fig.tight_layout()
    fig.savefig(f"work/fig1_a0{a0_slider.value:.2f}.png", dpi=300, bbox_inches='tight')
    fig
    return (fig,)


@app.cell
def _(mo, result):
    mo.md(
        f"""
        **均衡パラメータ**

        | 変数 | 値 |
        |------|-----|
        | $T_b$ | {result['Tb']:.4f} |
        | $T_p$ | {result['Tp']:.4f} |
        | $\\beta_{{b0}}$ | {result['betab']:.4f} |
        | $\\beta_{{p0}}$ | {result['betap']:.4f} |
        | $\\lambda_0$ | {result['lambda0']:.4f} |
        """
    )
    return


if __name__ == "__main__":
    app.run()
