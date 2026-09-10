import numpy as np
from scipy.linalg import expm

from .liouville import vec
from .pulses import (
    gaussian_spectrum,
    pulse_dressing_1,
    pulse_dressing_2,
    pulse_dressing_3,
)
from .resolvents import resolvent_action
from .response import (
    observable_bra,
    pathway_signs,
    finite_pulse_signal,
    impulsive_signal,
    short_pulse_rwa_signal,
    impulsive_rwa_signal,
)


def impulsive_spectrum(
    system,
    omega1,
    omega3,
    T,
    eta,
    rho0=None,
):
    """
    Evaluate the impulsive full-interaction third-order
    spectrum on a two-dimensional frequency grid.

    The omega1- and omega3-dependent pieces are evaluated
    separately and combined by matrix multiplication.
    """
    omega1 = np.asarray(
        omega1,
        dtype=float,
    )

    omega3 = np.asarray(
        omega3,
        dtype=float,
    )

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError(
            "omega1 and omega3 must be one-dimensional"
        )

    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    L = np.asarray(
        system.L,
        dtype=complex,
    )

    V = np.asarray(
        system.V,
        dtype=complex,
    )

    rho_vec = vec(rho0)

    mu_bra = observable_bra(
        system.dipole
    )

    n = L.shape[0]

    U_T = expm(
        L * T
    )

    # --------------------------------------------------------
    # omega1-dependent side
    #
    # R_i = V G1 V |rho0>>
    # --------------------------------------------------------

    right = np.empty(
        (n, omega1.size),
        dtype=complex,
    )

    for i, w1 in enumerate(omega1):

        state = V @ rho_vec

        state = resolvent_action(
            L_super=L,
            omega=w1,
            rhs=state,
            eta=eta,
        )

        state = V @ state

        right[:, i] = state

    # --------------------------------------------------------
    # omega3-dependent side
    #
    # L_j = <<mu| G3 V exp(LT)
    # --------------------------------------------------------

    left = np.empty(
        (omega3.size, n),
        dtype=complex,
    )

    for j, w3 in enumerate(omega3):

        bra_G3 = resolvent_action(
            L_super=L.T,
            omega=w3,
            rhs=mu_bra,
            eta=eta,
        )

        left[j, :] = (
            bra_G3
            @ V
            @ U_T
        )

    return left @ right

def response_grid(
    response_function,
    omega1,
    omega3,
    **kwargs,
):
    """
    Evaluate a frequency-domain response on a 2D spectral grid.

    Parameters
    ----------
    response_function : callable
        Function that accepts omega1 and omega3 as keyword arguments
        and returns a complex response value.
    omega1 : array_like
        Excitation-frequency grid.
    omega3 : array_like
        Detection-frequency grid.
    **kwargs
        Additional keyword arguments passed to response_function.

    Returns
    -------
    numpy.ndarray
        Complex response array with shape
        (len(omega3), len(omega1)).
    """
    omega1 = np.asarray(omega1, dtype=float)
    omega3 = np.asarray(omega3, dtype=float)

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError("omega1 and omega3 must be one-dimensional")

    spectrum = np.empty(
        (omega3.size, omega1.size),
        dtype=complex,
    )

    for j, w3 in enumerate(omega3):
        for i, w1 in enumerate(omega1):
            spectrum[j, i] = response_function(
                omega1=w1,
                omega3=w3,
                **kwargs,
            )

    return spectrum

def finite_pulse_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate a finite-pulse third-order spectrum on a
    two-dimensional frequency grid.

    The implementation caches all quantities that depend only
    on omega1 or omega3 and evaluates the final grid by matrix
    multiplication.

    Returns
    -------
    numpy.ndarray
        Complex spectrum with shape

            (len(omega3), len(omega1)).
    """
    omega1 = np.asarray(
        omega1,
        dtype=float,
    )

    omega3 = np.asarray(
        omega3,
        dtype=float,
    )

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError(
            "omega1 and omega3 must be one-dimensional"
        )

    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    L = np.asarray(
        system.L,
        dtype=complex,
    )

    V = np.asarray(
        system.V,
        dtype=complex,
    )

    rho_vec = vec(rho0)

    mu_bra = observable_bra(
        system.dipole
    )

    s1, s2, s3 = pathway_signs(
        pathway
    )

    n = L.shape[0]

    # Waiting-time propagation is identical everywhere
    # on the frequency grid.
    U_T = expm(
        L * T
    )

    # --------------------------------------------------------
    # Right-hand objects: depend only on omega1
    #
    # R_i = F2 V G1 V F1 |rho0>>
    # --------------------------------------------------------

    right = np.empty(
        (n, omega1.size),
        dtype=complex,
    )

    for i, w1 in enumerate(omega1):

        F1 = pulse_dressing_1(
            L_super=L,
            omega1=w1,
            omega_L=pulse1.omega_L,
            sigma=pulse1.sigma,
            E0=pulse1.E0,
            phase=pulse1.phase,
            sign=s1,
        )

        F2 = pulse_dressing_2(
            L_super=L,
            omega1=w1,
            omega_L=pulse2.omega_L,
            sigma=pulse2.sigma,
            E0=pulse2.E0,
            phase=pulse2.phase,
            sign=s2,
        )

        state = F1 @ rho_vec
        state = V @ state

        state = resolvent_action(
            L_super=L,
            omega=w1,
            rhs=state,
            eta=eta,
        )

        state = V @ state
        state = F2 @ state

        right[:, i] = state

    # --------------------------------------------------------
    # Left-hand objects: depend only on omega3
    #
    # L_j = <<mu| G3 V exp(LT) F3
    #
    # Instead of constructing G3 explicitly, solve the
    # transposed resolvent equation:
    #
    #     x = G3^T |mu>>
    #
    # so that
    #
    #     x^T = <<mu| G3.
    # --------------------------------------------------------

    left = np.empty(
        (omega3.size, n),
        dtype=complex,
    )

    for j, w3 in enumerate(omega3):

        F3 = pulse_dressing_3(
            L_super=L,
            omega3=w3,
            omega_L=pulse3.omega_L,
            sigma=pulse3.sigma,
            E0=pulse3.E0,
            phase=pulse3.phase,
            sign=s3,
        )

        bra_G3 = resolvent_action(
            L_super=L.T,
            omega=w3,
            rhs=mu_bra,
            eta=eta,
        )

        left[j, :] = (
            bra_G3
            @ V
            @ U_T
            @ F3
        )

    # --------------------------------------------------------
    # Entire 2D spectrum
    #
    # S[j, i] = L_j R_i
    # --------------------------------------------------------

    return left @ right


def short_pulse_rwa_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate the short-pulse RWA third-order spectrum.

    The molecular response is evaluated using the optimized
    impulsive RWA spectrum and then dressed by the three
    scalar Gaussian pulse spectra.
    """
    omega1 = np.asarray(
        omega1,
        dtype=float,
    )

    omega3 = np.asarray(
        omega3,
        dtype=float,
    )

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError(
            "omega1 and omega3 must be one-dimensional"
        )

    s1, s2, s3 = pathway_signs(
        pathway
    )

    molecular = impulsive_rwa_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway=pathway,
        rho0=rho0,
    )

    E1 = np.array(
        [
            gaussian_spectrum(
                omega=w1,
                omega_L=pulse1.omega_L,
                sigma=pulse1.sigma,
                E0=pulse1.E0,
                phase=pulse1.phase,
                sign=s1,
            )
            for w1 in omega1
        ],
        dtype=complex,
    )

    E2 = np.array(
        [
            gaussian_spectrum(
                omega=-w1,
                omega_L=pulse2.omega_L,
                sigma=pulse2.sigma,
                E0=pulse2.E0,
                phase=pulse2.phase,
                sign=s2,
            )
            for w1 in omega1
        ],
        dtype=complex,
    )

    E3 = np.array(
        [
            gaussian_spectrum(
                omega=w3,
                omega_L=pulse3.omega_L,
                sigma=pulse3.sigma,
                E0=pulse3.E0,
                phase=pulse3.phase,
                sign=s3,
            )
            for w3 in omega3
        ],
        dtype=complex,
    )

    return (
        E3[:, np.newaxis]
        * molecular
        * (E2 * E1)[np.newaxis, :]
    )

def impulsive_rwa_spectrum(
    system,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate the impulsive RWA third-order spectrum on a
    two-dimensional frequency grid.

    The omega1- and omega3-dependent parts are evaluated
    separately and combined by matrix multiplication.
    """
    omega1 = np.asarray(
        omega1,
        dtype=float,
    )

    omega3 = np.asarray(
        omega3,
        dtype=float,
    )

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError(
            "omega1 and omega3 must be one-dimensional"
        )

    if system.V_plus is None or system.V_minus is None:
        raise ValueError(
            "RWA calculations require dipole_plus and "
            "dipole_minus in SpectroscopySystem"
        )

    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    L = np.asarray(
        system.L,
        dtype=complex,
    )

    rho_vec = vec(rho0)

    mu_bra = observable_bra(
        system.dipole
    )

    s1, s2, s3 = pathway_signs(
        pathway
    )

    def interaction(sign):
        if sign == 1:
            return system.V_plus

        return system.V_minus

    V1 = interaction(s1)
    V2 = interaction(s2)
    V3 = interaction(s3)

    n = L.shape[0]

    U_T = expm(
        L * T
    )

    # --------------------------------------------------------
    # omega1-dependent right-hand side
    #
    # R_i = V2 G1 V1 |rho0>>
    # --------------------------------------------------------

    right = np.empty(
        (n, omega1.size),
        dtype=complex,
    )

    for i, w1 in enumerate(omega1):

        state = V1 @ rho_vec

        state = resolvent_action(
            L_super=L,
            omega=w1,
            rhs=state,
            eta=eta,
        )

        state = V2 @ state

        right[:, i] = state

    # --------------------------------------------------------
    # omega3-dependent left-hand side
    #
    # L_j = <<mu| G3 V3 exp(LT)
    # --------------------------------------------------------

    left = np.empty(
        (omega3.size, n),
        dtype=complex,
    )

    for j, w3 in enumerate(omega3):

        bra_G3 = resolvent_action(
            L_super=L.T,
            omega=w3,
            rhs=mu_bra,
            eta=eta,
        )

        left[j, :] = (
            bra_G3
            @ V3
            @ U_T
        )

    return left @ right
