import numpy as np

from .liouville import (
    liouvillian,
    interaction_superoperator,
    stationary_state,
)


class SpectroscopySystem:
    """
    Open quantum system for Liouville-space spectroscopy.

    Parameters
    ----------
    H : array_like
        System Hamiltonian.
    dipole : array_like
        Full dipole operator.
    collapse_ops : sequence of array_like, optional
        Lindblad collapse operators.
    rho0 : array_like, optional
        Initial density matrix.
    hbar : float, optional
        Reduced Planck constant. Default is 1.
    dipole_plus : array_like, optional
        Raising-frequency component of the dipole.
    dipole_minus : array_like, optional
        Lowering-frequency component of the dipole.
    """

    def __init__(
        self,
        H,
        dipole,
        collapse_ops=None,
        rho0=None,
        hbar=1.0,
        dipole_plus=None,
        dipole_minus=None,
    ):
        self.H = np.asarray(H, dtype=complex)
        self.dipole = np.asarray(dipole, dtype=complex)

        self.hbar = hbar
        self.d = self.H.shape[0]

        if self.H.shape != (self.d, self.d):
            raise ValueError("H must be a square matrix")

        if self.dipole.shape != (self.d, self.d):
            raise ValueError(
                "dipole must have the same dimensions as H"
            )

        if collapse_ops is None:
            self.collapse_ops = []
        else:
            self.collapse_ops = [
                np.asarray(C, dtype=complex)
                for C in collapse_ops
            ]

        if rho0 is None:
            self.rho0 = None
        else:
            self.rho0 = np.asarray(rho0, dtype=complex)

        self.L = liouvillian(
            self.H,
            collapse_ops=self.collapse_ops,
            hbar=self.hbar,
        )

        self.V = interaction_superoperator(
            self.dipole,
            hbar=self.hbar,
        )

        self.dipole_plus = None
        self.dipole_minus = None

        self.V_plus = None
        self.V_minus = None

        if dipole_plus is not None:
            self.dipole_plus = np.asarray(
                dipole_plus,
                dtype=complex,
            )

            self.V_plus = interaction_superoperator(
                self.dipole_plus,
                hbar=self.hbar,
            )

        if dipole_minus is not None:
            self.dipole_minus = np.asarray(
                dipole_minus,
                dtype=complex,
            )

            self.V_minus = interaction_superoperator(
                self.dipole_minus,
                hbar=self.hbar,
            )

    def stationary_state(self):
        """
        Return the stationary state of the field-free Liouvillian.
        """
        return stationary_state(
            self.L,
            self.d,
        )
