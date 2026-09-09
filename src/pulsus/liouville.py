import numpy as np


def vec(rho):
    """
    Column-stacking vectorization of a density matrix.

    Parameters
    ----------
    rho : array_like
        Square matrix to vectorize.

    Returns
    -------
    numpy.ndarray
        Column-stacked vector representation.
    """
    rho = np.asarray(rho, dtype=complex)
    return rho.reshape(-1, order="F")


def unvec(rho_vec, d):
    """
    Inverse of column-stacking vectorization.

    Parameters
    ----------
    rho_vec : array_like
        Vectorized density matrix.
    d : int
        Hilbert-space dimension.

    Returns
    -------
    numpy.ndarray
        Density matrix with shape (d, d).
    """
    rho_vec = np.asarray(rho_vec, dtype=complex)
    return rho_vec.reshape((d, d), order="F")


def hamiltonian_liouvillian(H, hbar=1.0):
    """
    Construct the Hamiltonian part of the Liouvillian.

    The convention is

        L_H |rho>> = -(i / hbar) [H, rho]

    using column-stacking vectorization.

    Parameters
    ----------
    H : array_like
        System Hamiltonian.
    hbar : float, optional
        Reduced Planck constant. Default is 1.

    Returns
    -------
    numpy.ndarray
        Hamiltonian Liouvillian superoperator.
    """
    H = np.asarray(H, dtype=complex)

    d = H.shape[0]
    I = np.eye(d, dtype=complex)

    return (
        -1j / hbar
        * (
            np.kron(I, H)
            - np.kron(H.T, I)
        )
    )


def interaction_superoperator(mu, hbar=1.0):
    """
    Construct the light-matter interaction superoperator.

    The convention is

        V |rho>> = (i / hbar) [mu, rho]

    using column-stacking vectorization.

    Parameters
    ----------
    mu : array_like
        Interaction operator, typically the dipole operator.
    hbar : float, optional
        Reduced Planck constant. Default is 1.

    Returns
    -------
    numpy.ndarray
        Interaction superoperator.
    """
    mu = np.asarray(mu, dtype=complex)

    d = mu.shape[0]
    I = np.eye(d, dtype=complex)

    return (
        1j / hbar
        * (
            np.kron(I, mu)
            - np.kron(mu.T, I)
        )
    )


def dissipator_super(L, I=None):
    """
    Liouville-space Lindblad dissipator using column stacking.

    The dissipator corresponds to

        D[L](rho) = L rho L^dagger
                    - 1/2 {L^dagger L, rho}.

    Parameters
    ----------
    L : array_like
        Lindblad jump operator.
    I : array_like, optional
        Identity matrix. If omitted, it is constructed automatically.

    Returns
    -------
    numpy.ndarray
        Liouville-space dissipator superoperator.
    """
    L = np.asarray(L, dtype=complex)

    if I is None:
        I = np.eye(L.shape[0], dtype=complex)
    else:
        I = np.asarray(I, dtype=complex)

    LdagL = L.conj().T @ L

    return (
        np.kron(L.conj(), L)
        - 0.5 * np.kron(I, LdagL)
        - 0.5 * np.kron(LdagL.T, I)
    )


def liouvillian(H, collapse_ops=None, hbar=1.0):
    """
    Construct the full field-free Liouvillian.

    Parameters
    ----------
    H : array_like
        System Hamiltonian.
    collapse_ops : sequence of array_like, optional
        Lindblad collapse operators. Rates should already be
        included in the operators.
    hbar : float, optional
        Reduced Planck constant. Default is 1.

    Returns
    -------
    numpy.ndarray
        Full Liouvillian superoperator.
    """
    L_super = hamiltonian_liouvillian(H, hbar=hbar)

    if collapse_ops is None:
        return L_super

    for C in collapse_ops:
        L_super = L_super + dissipator_super(C)

    return L_super


def stationary_state(L_super, d):
    """
    Return the normalized stationary density matrix associated
    with the Liouvillian eigenvalue closest to zero.

    Parameters
    ----------
    L_super : array_like
        Liouvillian superoperator.
    d : int
        Hilbert-space dimension.

    Returns
    -------
    numpy.ndarray
        Normalized stationary density matrix.
    """
    L_super = np.asarray(L_super, dtype=complex)

    eigenvalues, eigenvectors = np.linalg.eig(L_super)
    zero_index = np.argmin(np.abs(eigenvalues))

    rho_ss_vec = eigenvectors[:, zero_index]
    rho_ss = unvec(rho_ss_vec, d)

    rho_ss = rho_ss / np.trace(rho_ss)

    # Remove small numerical non-Hermitian components.
    rho_ss = 0.5 * (rho_ss + rho_ss.conj().T)

    rho_ss = rho_ss / np.trace(rho_ss)

    return rho_ss
