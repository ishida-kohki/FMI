"""Filament Merging Instability (FMI) analysis package.

Based on Vanthieghem et al., Phys. Plasmas 25, 072115 (2018).
"""

from fmi import (
    equilibrium,
    floquet,
    mach_parameters,
    snr_equilibrium,
    snr_linear,
    snr_linear_relative,
    wi_saturated_temperature,
)

__all__ = [
    "equilibrium",
    "floquet",
    "mach_parameters",
    "snr_equilibrium",
    "snr_linear",
    "snr_linear_relative",
    "wi_saturated_temperature",
]
