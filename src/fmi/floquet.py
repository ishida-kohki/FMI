"""Floquet analysis for the FMI dispersion relation.

Vanthieghem et al., Phys. Plasmas 25, 072115 (2018), Appendix A and Sec. III A.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp


def get_N_matrix(y, omega, kx, params):
    c, q, eps, beta0, gamma0 = params['c'], params['q'], params['eps'], params['beta0'], params['gamma0']
    d0, p0, G_ad, mass = params['d0'], params['p0'], params['Gamma_ad'], params['m']

    B0z_const, E0y, f_B, f_E, f_a = params['bg_func'](y)

    C_ey_bz = kx * c / omega
    C_ey_Vy = np.zeros(4, dtype=np.complex128)
    for b in range(4):
        C_ey_Vy[b] = -1j * 4.0 * np.pi * c / (omega * B0z_const) * q[b] * beta0[b] * d0[b] * f_a[b]

    W = np.zeros(4, dtype=np.complex128)
    C_Vx_ex = np.zeros(4, dtype=np.complex128)
    C_Vx_P  = np.zeros(4, dtype=np.complex128)
    C_Vx_Vy = np.zeros(4, dtype=np.complex128)
    K_a  = np.zeros(4, dtype=np.complex128)
    C_D_ex = np.zeros(4, dtype=np.complex128)
    C_D_P  = np.zeros(4, dtype=np.complex128)
    C_D_Vy = np.zeros(4, dtype=np.complex128)

    for a in range(4):
        mass_term = mass * c**2 * d0[a] / gamma0[a]
        h0_a = p0[a] * G_ad[a] / (G_ad[a] - 1.0) + mass_term

        W[a] = (gamma0[a]**2 / c**2) * h0_a * f_a[a] * beta0[a] * (-1j * omega + 1j * eps[a] * beta0[a] * c * kx)

        C_Vx_P[a]  = (-1j * kx * p0[a] + 1j * eps[a] * omega * p0[a] * beta0[a] / c) / W[a]
        C_Vx_ex[a] = (q[a] * d0[a] * B0z_const * f_a[a] * (1.0 - beta0[a]**2)) / W[a]
        C_Vx_Vy[a] = (q[a] * d0[a] * beta0[a] * f_a[a] * (B0z_const * f_B - eps[a] * beta0[a] * E0y * f_E)) / W[a]

        K_a[a]  = eps[a] * f_a[a] * (beta0[a]**2) / (1.0 - beta0[a]**2)
        C_D_ex[a] = K_a[a] * C_Vx_ex[a]
        C_D_P[a]  = 1.0 / G_ad[a] + K_a[a] * C_Vx_P[a]
        C_D_Vy[a] = K_a[a] * C_Vx_Vy[a]

    N = np.zeros((10, 10), dtype=np.complex128)

    # eq(A1)
    N[0, 1] = 1j * c / omega * (kx**2 - (omega / c)**2)
    for b in range(4):
        N[0, 6 + b] = 1j * kx * C_ey_Vy[b]

    # eq(A2)
    N[1, 0] = -1j * omega / c
    for a in range(4):
        S_a = (4.0 * np.pi / B0z_const) * q[a] * beta0[a] * d0[a]
        N[1, 0]     += S_a * (f_a[a] * C_Vx_ex[a] + eps[a] * C_D_ex[a])
        N[1, 2 + a]  = S_a * (f_a[a] * C_Vx_P[a]  + eps[a] * C_D_P[a])
        N[1, 6 + a]  = S_a * (f_a[a] * C_Vx_Vy[a] + eps[a] * C_D_Vy[a])

    # eq(A3)
    for a in range(4):
        H_a = (q[a] * d0[a] * B0z_const) / p0[a]
        Q_a = -eps[a] * beta0[a] * f_B + (E0y / B0z_const) * f_E

        mass_term = mass * c**2 * d0[a] / gamma0[a]
        h0_a = p0[a] * G_ad[a] / (G_ad[a] - 1.0) + mass_term
        U_a = (gamma0[a]**2 * beta0[a] / (c * p0[a])) * h0_a * f_a[a]

        N[2 + a, 0]     = H_a * (Q_a * C_D_ex[a] - beta0[a] * f_a[a] * f_B * C_Vx_ex[a])
        N[2 + a, 1]     = H_a * (f_a[a] * C_ey_bz - eps[a] * beta0[a] * f_a[a])
        N[2 + a, 2 + a] = H_a * (Q_a * C_D_P[a]  - beta0[a] * f_a[a] * f_B * C_Vx_P[a])

        for b in range(4):
            if b != a:
                N[2 + a, 6 + b] = H_a * f_a[a] * C_ey_Vy[b]

        N[2 + a, 6 + a] = H_a * (Q_a * C_D_Vy[a] + f_a[a] * C_ey_Vy[a] - beta0[a] * f_a[a] * f_B * C_Vx_Vy[a]) \
                          - U_a * (-1j * omega + 1j * kx * c * eps[a] * beta0[a])

    # eq(A4)
    for a in range(4):
        R_a    = 1j * omega / (c * beta0[a] * f_a[a]) - 1j * kx * eps[a] / f_a[a]
        beta0x = eps[a] * beta0[a]
        dln_fa = (q[a] * d0[a] / p0[a]) * (E0y * f_E - beta0x * B0z_const * f_B)

        N[6 + a, 0]     = R_a * C_D_ex[a] - 1j * kx * C_Vx_ex[a]
        N[6 + a, 2 + a] = R_a * C_D_P[a]  - 1j * kx * C_Vx_P[a]
        N[6 + a, 6 + a] = R_a * C_D_Vy[a] - 1j * kx * C_Vx_Vy[a] - dln_fa

    return N


def evaluate_and_verify_floquet(omega, kx, params):
    """Compute Floquet multipliers and return Bloch wavenumbers k_y.

    Integrates dX/dy = N(y) X forward from lambda0/2 to lambda0 and backward
    from lambda0/2 to 0, then forms B = inv(X_0) @ X_L.  Eigenvalues with
    |rho| ~ 1 are real (propagating) Bloch modes; their wavenumber is
    k_y = arg(rho) / lambda0.

    Returns
    -------
    list of float
        k_y values (in physical units, not normalized).
    """
    lambda0 = params['lambda0']
    X_init_r = np.concatenate([np.eye(10, dtype=np.complex128).real.flatten(),
                                np.eye(10, dtype=np.complex128).imag.flatten()])

    def ode_system_real(y, X_flat_real):
        X_R = X_flat_real[:100].reshape((10, 10))
        X_I = X_flat_real[100:].reshape((10, 10))
        N = get_N_matrix(y, omega, kx, params)
        dX_R = N.real @ X_R - N.imag @ X_I
        dX_I = N.real @ X_I + N.imag @ X_R
        return np.concatenate([dX_R.flatten(), dX_I.flatten()])

    sol_fwd = solve_ivp(ode_system_real, [lambda0 / 2.0, lambda0],  X_init_r, method='Radau', rtol=1e-8, atol=1e-8)
    sol_bwd = solve_ivp(ode_system_real, [lambda0 / 2.0, 0.0],     X_init_r, method='Radau', rtol=1e-8, atol=1e-8)

    X_L = sol_fwd.y[:100, -1].reshape((10, 10)) + 1j * sol_fwd.y[100:, -1].reshape((10, 10))
    X_0 = sol_bwd.y[:100, -1].reshape((10, 10)) + 1j * sol_bwd.y[100:, -1].reshape((10, 10))

    B = la.inv(X_0) @ X_L
    eigvals = la.eigvals(B)

    unfolded_kys = []
    for rho in eigvals:
        if abs(abs(rho) - 1.0) < 1e-4:
            theta = np.angle(rho)
            ky = theta / lambda0
            unfolded_kys.append(ky)

    return unfolded_kys


def scan_dispersion(kx, params, Omega_p, k0_paper, gamma_list):
    """Scan growth rate Gamma and return dispersion relation with Brillouin zone unfolding.

    Wraps ``evaluate_and_verify_floquet`` over ``gamma_list`` and applies the
    unfolding logic: when the branch hits the gap at k0/2 and starts folding back,
    the returning branch is mirrored as k0_paper - k to restore the full curve.

    Parameters
    ----------
    kx : float
        Longitudinal wavenumber (physical units).
    params : dict
        Equilibrium parameter dictionary.
    Omega_p : float
        Plasma frequency.
    k0_paper : float
        Fundamental wavenumber k0 normalized by Omega_p/c.
    gamma_list : array-like
        Growth rates Gamma to scan (decreasing order recommended).

    Returns
    -------
    ky_plot : list of float
        k_y values normalized by Omega_p/c.
    gamma_plot : list of float
        Corresponding growth rates Gamma in units of Omega_p.
    """
    ky_plot = []
    gamma_plot = []

    has_hit_gap = False
    prev_max_ky = 0.0

    for gamma_paper in gamma_list:
        omega_test = gamma_paper * Omega_p * 1j
        unfolded_kys = evaluate_and_verify_floquet(omega_test, kx, params)

        unique_kys = []
        for k in unfolded_kys:
            k_norm = k / Omega_p
            if not any(abs(k_norm - uk) < 1e-4 for uk in unique_kys):
                unique_kys.append(k_norm)

        if not unique_kys:
            continue

        current_max_ky = max(unique_kys)

        # ギャップ（k0/2）に衝突し、位相の巻き戻し（値の減少）が始まった瞬間を検知
        if not has_hit_gap and prev_max_ky > 0 and current_max_ky < prev_max_ky - 1e-4:
            has_hit_gap = True

        for k in unique_kys:
            # 折り返しフラグが立っており、かつ k が「戻ってきた枝（大きい方の値）」なら対称移動
            if has_hit_gap and abs(k - current_max_ky) < 1e-4:
                ky_unfolded = k0_paper - k
                ky_plot.append(ky_unfolded)
                gamma_plot.append(gamma_paper)
            else:
                # それ以外（左側の主枝）はそのままプロット
                ky_plot.append(k)
                gamma_plot.append(gamma_paper)

        prev_max_ky = current_max_ky

    return ky_plot, gamma_plot
