import numpy as np


def vec(rho):
    """
    Column-stacking vectorization.
    """
    return rho.reshape(-1, order="F")


def unvec(rho_vec, d):
    """
    Inverse of column-stacking vectorization.
    """
    return rho_vec.reshape((d, d), order="F")


def dissipator_super(L, I):
    """
    Liouville-space Lindblad dissipator.
    """

    LdagL = L.conj().T @ L

    return (
        np.kron(L.conj(), L)
        - 0.5 * np.kron(I, LdagL)
        - 0.5 * np.kron(LdagL.T, I)
    )


def build_dimer(
    omega_a=1.0,
    omega_b=1.2,
    J=0.10,
    gamma_a=0.05,
    gamma_b=0.05,
    nbar_a=0.05,
    nbar_b=0.05,
    mu_a=1.0,
    mu_b=1.0,
    hbar=1.0,
):

    # --------------------------------------------------------
    # Single-TLS basis and operators
    # --------------------------------------------------------

    g = np.array([1.0, 0.0], dtype=complex)
    e = np.array([0.0, 1.0], dtype=complex)

    I2 = np.eye(2, dtype=complex)

    sigma_plus = np.outer(e, g.conj())
    sigma_minus = np.outer(g, e.conj())

    # --------------------------------------------------------
    # Composite basis
    #
    # Tensor convention:
    # H = H_b ⊗ H_a
    #
    # Numerical basis:
    # |gg>, |eg>, |ge>, |ee>
    # --------------------------------------------------------

    gg = np.kron(g, g)
    eg = np.kron(g, e)
    ge = np.kron(e, g)
    ee = np.kron(e, e)

    sigma_a_plus = np.kron(I2, sigma_plus)
    sigma_a_minus = np.kron(I2, sigma_minus)

    sigma_b_plus = np.kron(sigma_plus, I2)
    sigma_b_minus = np.kron(sigma_minus, I2)

    # --------------------------------------------------------
    # Hamiltonian
    # --------------------------------------------------------

    n_a = sigma_a_plus @ sigma_a_minus
    n_b = sigma_b_plus @ sigma_b_minus

    exchange = (
        sigma_a_plus @ sigma_b_minus
        + sigma_a_minus @ sigma_b_plus
    )

    H_S = (
        hbar * omega_a * n_a
        + hbar * omega_b * n_b
        + hbar * J * exchange
    )

    # --------------------------------------------------------
    # Dipole
    # --------------------------------------------------------

    mu = (
        mu_a * (sigma_a_plus + sigma_a_minus)
        + mu_b * (sigma_b_plus + sigma_b_minus)
    )

    # --------------------------------------------------------
    # Thermal jump operators
    # --------------------------------------------------------

    L_a_down = (
        np.sqrt(gamma_a * (nbar_a + 1.0))
        * sigma_a_minus
    )

    L_a_up = (
        np.sqrt(gamma_a * nbar_a)
        * sigma_a_plus
    )

    L_b_down = (
        np.sqrt(gamma_b * (nbar_b + 1.0))
        * sigma_b_minus
    )

    L_b_up = (
        np.sqrt(gamma_b * nbar_b)
        * sigma_b_plus
    )

    # --------------------------------------------------------
    # Liouvillian
    # --------------------------------------------------------

    d = H_S.shape[0]
    I4 = np.eye(d, dtype=complex)

    L_H = (
        -1j / hbar
        * (
            np.kron(I4, H_S)
            - np.kron(H_S.T, I4)
        )
    )

    D_a_down = dissipator_super(L_a_down, I4)
    D_a_up = dissipator_super(L_a_up, I4)

    D_b_down = dissipator_super(L_b_down, I4)
    D_b_up = dissipator_super(L_b_up, I4)

    L_super = (
        L_H
        + D_a_down
        + D_a_up
        + D_b_down
        + D_b_up
    )

    # --------------------------------------------------------
    # Light-matter interaction superoperator
    # --------------------------------------------------------

    V = (
        1j / hbar
        * (
            np.kron(I4, mu)
            - np.kron(mu.T, I4)
        )
    )

    # --------------------------------------------------------
    # Initial state
    # --------------------------------------------------------

    rho0 = np.outer(gg, gg.conj())

    # --------------------------------------------------------
    # Return model
    # --------------------------------------------------------

    return {
        "hbar": hbar,
        "d": d,

        "omega_a": omega_a,
        "omega_b": omega_b,
        "J": J,

        "gamma_a": gamma_a,
        "gamma_b": gamma_b,
        "nbar_a": nbar_a,
        "nbar_b": nbar_b,

        "gg": gg,
        "eg": eg,
        "ge": ge,
        "ee": ee,

        "sigma_a_plus": sigma_a_plus,
        "sigma_a_minus": sigma_a_minus,
        "sigma_b_plus": sigma_b_plus,
        "sigma_b_minus": sigma_b_minus,

        "H_S": H_S,
        "mu": mu,

        "L_a_down": L_a_down,
        "L_a_up": L_a_up,
        "L_b_down": L_b_down,
        "L_b_up": L_b_up,

        "L_super": L_super,
        "V": V,

        "rho0": rho0,
        "I4": I4,
    }

def stationary_state(L_super, d):
    """
    Return the normalized stationary density matrix associated
    with the Liouvillian eigenvalue closest to zero.
    """

    eigenvalues, eigenvectors = np.linalg.eig(L_super)

    zero_index = np.argmin(np.abs(eigenvalues))

    rho_ss_vec = eigenvectors[:, zero_index]
    rho_ss = unvec(rho_ss_vec, d)

    # Normalize to unit trace
    rho_ss = rho_ss / np.trace(rho_ss)

    # Remove tiny numerical non-Hermitian components
    rho_ss = 0.5 * (rho_ss + rho_ss.conj().T)

    # Renormalize after symmetrization
    rho_ss = rho_ss / np.trace(rho_ss)

    return rho_ss
