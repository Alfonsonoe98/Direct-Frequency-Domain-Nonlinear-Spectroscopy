"""Shared model for the repository tutorials; not part of the pulsus API."""

import numpy as np
import pulsus


def build_dimer():
    """Return the tutorial dimer and its stationary initial state (hbar=1)."""
    identity = np.eye(2, dtype=complex)
    lower = np.array([[0, 1], [0, 0]], dtype=complex)
    a = np.kron(identity, lower)
    b = np.kron(lower, identity)
    H = (a.conj().T @ a + 1.2 * (b.conj().T @ b)
         + 0.1 * (a.conj().T @ b + b.conj().T @ a))
    mu_plus = a.conj().T + b.conj().T
    gamma, nbar = 0.05, 0.05
    collapse_ops = [
        np.sqrt(gamma * (nbar + 1)) * a,
        np.sqrt(gamma * nbar) * a.conj().T,
        np.sqrt(gamma * (nbar + 1)) * b,
        np.sqrt(gamma * nbar) * b.conj().T,
    ]
    system = pulsus.SpectroscopySystem(
        H=H, dipole=mu_plus + mu_plus.conj().T,
        collapse_ops=collapse_ops, hbar=1.0,
        dipole_plus=mu_plus, dipole_minus=mu_plus.conj().T,
    )
    return system, system.stationary_state()
