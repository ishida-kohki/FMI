"""Equilibrium profiles for the FMI (Filament Merging Instability).

Implements the analytical equilibrium of Vanthieghem et al.,
Phys. Plasmas 25, 072115 (2018).
"""

from __future__ import annotations

import numpy as np
import scipy.special as sp


def build_equilibrium(xi: float = 5.0, T0: float = 1.0, gamma0: float = 10.0):
    """Build equilibrium parameters for a pair-plasma filament lattice.

    Parameters
    ----------
    xi:
        Magnetization parameter.
    T0:
        Reference temperature (normalized).
    gamma0:
        Lorentz factor of the beam.

    Returns
    -------
    params:
        Dictionary of equilibrium parameters consumed by the Floquet solver.
    Omega_p:
        Plasma frequency (= exp(xi/2) in normalized units).
    k0_paper:
        Fundamental wavenumber k_0 normalized by Omega_p/c.
    """
    beta0 = np.sqrt(1.0 - 1.0 / gamma0**2)
    m1 = np.tanh(xi / 2.0) ** 2
    K_m1 = sp.ellipk(m1)
    u1 = 2.0 * beta0 * gamma0 / np.sqrt(T0) * np.cosh(xi / 2.0)
    u = u1 * np.tanh(xi / 2.0)
    lambda0 = 4.0 * K_m1 / u1
    B0z_const = 2.0 * u * T0 / (gamma0 * beta0)

    # Species ordering (paper Appendix A):
    # Fluid index a=0,1,2,3 → eps=(+1,-1,-1,+1), q=(-1,-1,+1,+1)
    # i.e. right e-, left e-, left e+, right e+
    eps_arr = np.array([1.0, -1.0, -1.0, 1.0])
    q_arr = np.array([-1.0, -1.0, 1.0, 1.0])

    def bg_func(y):
        """Return background fields at position y (scalar or array)."""
        phase = u1 * (y - lambda0 / 4.0)
        sn, cn, dn, _ = sp.ellipj(phase, m1)
        f_B = cn / dn  # B0z(y) / B0z_const
        tanh_xi2 = np.tanh(xi / 2.0)
        n = len(y) if isinstance(y, np.ndarray) else 1
        fa = np.zeros((4, n), dtype=np.float64)
        for a in range(4):
            fa[a] = (dn / (1.0 + eps_arr[a] * q_arr[a] * tanh_xi2 * sn)) ** 2
        if not isinstance(y, np.ndarray):
            fa = fa[:, 0]
            f_B = float(f_B)
        return B0z_const, 0.0, f_B, 0.0, fa

    n0_val = 1.0 / (4.0 * np.pi)
    d0_val = gamma0 * n0_val
    p0_val = n0_val * T0

    params = {
        "c": 1.0,
        "m": 1.0,
        "q": q_arr,
        "eps": eps_arr,
        "beta0": np.full(4, beta0),
        "gamma0": np.full(4, gamma0),
        "d0": np.full(4, d0_val),
        "p0": np.full(4, p0_val),
        "Gamma_ad": np.full(4, 4.0 / 3.0),
        "lambda0": lambda0,
        "bg_func": bg_func,
    }

    Omega_p = np.exp(xi / 2.0)
    k0_paper = (2.0 * np.pi / lambda0) / Omega_p
    return params, Omega_p, k0_paper
