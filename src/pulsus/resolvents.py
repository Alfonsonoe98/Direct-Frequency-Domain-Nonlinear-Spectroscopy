import numpy as np
from scipy.linalg import solve


def resolvent_action(L_super, omega, rhs, eta):
    """
    Apply the frequency-domain Liouvillian resolvent to a vector
    or matrix of right-hand sides.

    The resolvent is

        G_eta(omega) = [(eta - i*omega) I - L]^{-1}.

    Rather than constructing the inverse explicitly, this function
    solves

        [(eta - i*omega) I - L] x = rhs.

    Parameters
    ----------
    L_super : array_like
        Liouvillian superoperator.
    omega : float
        Frequency at which the resolvent is evaluated.
    rhs : array_like
        Right-hand side vector or matrix.
    eta : float
        Positive convergence/broadening parameter.

    Returns
    -------
    numpy.ndarray
        Result of G_eta(omega) acting on rhs.
    """
    L_super = np.asarray(L_super, dtype=complex)
    rhs = np.asarray(rhs, dtype=complex)

    n = L_super.shape[0]
    I = np.eye(n, dtype=complex)

    A = (eta - 1j * omega) * I - L_super

    return solve(A, rhs)


def resolvent_matrix(L_super, omega, eta):
    """
    Construct the full frequency-domain Liouvillian resolvent matrix.

    Parameters
    ----------
    L_super : array_like
        Liouvillian superoperator.
    omega : float
        Frequency at which the resolvent is evaluated.
    eta : float
        Positive convergence/broadening parameter.

    Returns
    -------
    numpy.ndarray
        Full resolvent matrix G_eta(omega).
    """
    L_super = np.asarray(L_super, dtype=complex)

    n = L_super.shape[0]
    I = np.eye(n, dtype=complex)

    return resolvent_action(
        L_super,
        omega=omega,
        rhs=I,
        eta=eta,
    )
