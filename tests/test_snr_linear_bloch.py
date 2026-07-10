"""Fourier-Bloch 線形固有値ソルバの回帰テスト。"""

import numpy as np
import pytest

from fmi.snr_linear import BackgroundProfiles
from fmi.snr_linear_bloch import build_operator_bloch


def _uniform_single_electron(M: int = 24, L: float = 2 * np.pi) -> BackgroundProfiles:
    y = np.linspace(0.0, L, M, endpoint=False)
    return BackgroundProfiles(
        y=y, n={"e": np.ones(M)}, beta={"e": 0.0}, T={"e": 0.0},
        mass={"e": 1.0}, charge={"e": -1.0}, B0z=np.zeros(M),
        lambda0=L, n_periods=1,
    )


def test_bloch_lightwave_K0():
    """K=0 で Bloch 演算子は光波/Langmuir 分散 ω²=1+kx² を再現。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5, 1.3):
        w2 = np.linalg.eigvals(build_operator_bloch(bg, kx, K=0.0)) ** 2
        assert np.min(np.abs(w2 - (1.0 + kx**2))) < 1e-8
        assert np.min(np.abs(w2 - 1.0)) < 1e-8


def test_bloch_shift_adds_ky():
    """K≠0 は面内波数 ky=K として入り ω²=1+kx²+K² を与える（iK シフトの検証）。"""
    bg = _uniform_single_electron()
    for kx in (0.0, 0.5):
        for K in (0.3, 0.7):
            w2 = np.linalg.eigvals(build_operator_bloch(bg, kx, K=K)) ** 2
            assert np.min(np.abs(w2 - (1.0 + kx**2 + K**2))) < 1e-8
