import numpy as np
from scipy.linalg import expm

from .liouville import vec
from .resolvents import resolvent_action

from .pulses import (
    pulse_dressing_1,
    pulse_dressing_2,
    pulse_dressing_3,
    gaussian_spectrum,
)

def short_pulse_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma1,
    sigma2,
    sigma3,
    pathway="NR",
    E01=1.0,
    E02=1.0,
    E03=1.0,
    phase1=0.0,
    phase2=0.0,
    phase3=0.0,
):
    """
    Evaluate the short-pulse RWA response.

    Molecular evolution during each pulse is neglected, while
    the Gaussian spectral amplitudes are retained.
    """
    s1, s2, s3 = pathway_signs(pathway)

    molecular = impulsive_rwa_response(
        L_super=L_super,
        V_plus=V_plus,
        V_minus=V_minus,
        observable=observable,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway=pathway,
    )

    E1 = gaussian_spectrum(
        omega=omega1,
        omega_L=omega_L1,
        sigma=sigma1,
        E0=E01,
        phase=phase1,
        sign=s1,
    )

    # Pulse 2 carries the conjugate excitation-frequency argument.
    E2 = gaussian_spectrum(
        omega=-omega1,
        omega_L=omega_L2,
        sigma=sigma2,
        E0=E02,
        phase=phase2,
        sign=s2,
    )

    E3 = gaussian_spectrum(
        omega=omega3,
        omega_L=omega_L3,
        sigma=sigma3,
        E0=E03,
        phase=phase3,
        sign=s3,
    )

    return E3 * E2 * E1 * molecular

def impulsive_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
):
    """
    Evaluate the impulsive third-order response in the
    rotating-wave approximation.
    """
    s1, s2, s3 = pathway_signs(pathway)

    interactions = {
        1: np.asarray(V_plus, dtype=complex),
        -1: np.asarray(V_minus, dtype=complex),
    }

    V1 = interactions[s1]
    V2 = interactions[s2]
    V3 = interactions[s3]

    rho0 = np.asarray(rho0, dtype=complex)

    if rho0.ndim == 2:
        state = vec(rho0)
    elif rho0.ndim == 1:
        state = rho0.copy()
    else:
        raise ValueError(
            "rho0 must be a density matrix or vectorized state"
        )

    state = V1 @ state

    state = resolvent_action(
        L_super,
        omega=omega1,
        rhs=state,
        eta=eta,
    )

    state = V2 @ state

    state = expm(L_super * T) @ state

    state = V3 @ state

    state = resolvent_action(
        L_super,
        omega=omega3,
        rhs=state,
        eta=eta,
    )

    bra = observable_bra(observable)

    return bra @ state

def observable_bra(observable):
    """
    Construct the Liouville-space bra corresponding to an observable.

    For column-stacking vectorization,

        <A>> |rho>> = Tr(A rho).

    Parameters
    ----------
    observable : array_like
        Detection observable.

    Returns
    -------
    numpy.ndarray
        Liouville-space bra represented as a one-dimensional array.
    """
    observable = np.asarray(observable, dtype=complex)

    return vec(observable.conj().T).conj()


def third_order_response(
    L_super,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    F1,
    F2,
    F3,
):
    """
    Evaluate the finite-pulse third-order frequency-domain response

        P^(3)(omega3, T, omega1)
        =
        <<mu|
        G(omega3) V exp(L T)
        F3 F2 V G(omega1) V F1
        |rho0>>.

    Parameters
    ----------
    L_super : array_like
        Field-free Liouvillian.
    V : array_like
        Light-matter interaction superoperator.
    observable : array_like
        Detection observable, typically the dipole operator.
    rho0 : array_like
        Initial density matrix or vectorized density matrix.
    omega1, omega3 : float
        Excitation and detection frequencies.
    T : float
        Waiting time.
    eta : float
        Resolvent convergence/broadening parameter.
    F1, F2, F3 : array_like
        Pulse-dressing operators.

    Returns
    -------
    complex
        Third-order response.
    """
    L_super = np.asarray(L_super, dtype=complex)
    V = np.asarray(V, dtype=complex)

    F1 = np.asarray(F1, dtype=complex)
    F2 = np.asarray(F2, dtype=complex)
    F3 = np.asarray(F3, dtype=complex)

    rho0 = np.asarray(rho0, dtype=complex)

    if rho0.ndim == 2:
        state = vec(rho0)
    elif rho0.ndim == 1:
        state = rho0.copy()
    else:
        raise ValueError("rho0 must be a density matrix or vectorized state")

    # First pulse dressing
    state = F1 @ state

    # First interaction
    state = V @ state

    # First frequency-domain interval
    state = resolvent_action(
        L_super,
        omega=omega1,
        rhs=state,
        eta=eta,
    )

    # Second interaction
    state = V @ state

    # Finite-pulse dressing associated with pulses 2 and 3
    state = F2 @ state
    state = F3 @ state

    # Waiting-time evolution
    state = expm(L_super * T) @ state

    # Third interaction
    state = V @ state

    # Detection-frequency interval
    state = resolvent_action(
        L_super,
        omega=omega3,
        rhs=state,
        eta=eta,
    )

    # Detection
    bra = observable_bra(observable)

    return bra @ state

def pathway_signs(pathway):
    """
    Return the signed-frequency components for a standard
    third-order rephasing or nonrephasing signal.

    Parameters
    ----------
    pathway : {"NR", "R"}
        Nonrephasing ("NR") or rephasing ("R").

    Returns
    -------
    tuple
        (s1, s2, s3)
    """
    pathway = pathway.upper()

    if pathway == "NR":
        return 1, -1, 1

    if pathway == "R":
        return -1, 1, 1

    raise ValueError("pathway must be 'NR' or 'R'")


def finite_pulse_response(
    L_super,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma1,
    sigma2,
    sigma3,
    pathway="NR",
    E01=1.0,
    E02=1.0,
    E03=1.0,
    phase1=0.0,
    phase2=0.0,
    phase3=0.0,
):
    """
    Evaluate the finite-pulse third-order response for a
    rephasing or nonrephasing field-sign sector.

    The excitation frequency omega1 is signed:

        omega1 > 0 for NR
        omega1 < 0 for R

    Parameters
    ----------
    L_super : array_like
        Field-free Liouvillian.
    V : array_like
        Full light-matter interaction superoperator.
    observable : array_like
        Detection observable, typically the dipole operator.
    rho0 : array_like
        Initial density matrix or vectorized state.
    omega1 : float
        Signed excitation frequency.
    omega3 : float
        Detection frequency.
    T : float
        Waiting time.
    eta : float
        Resolvent broadening parameter.
    omega_L1, omega_L2, omega_L3 : float
        Carrier frequencies of the three pulses.
    sigma1, sigma2, sigma3 : float
        Gaussian temporal widths of the three pulses.
    pathway : {"NR", "R"}, optional
        Nonrephasing or rephasing field-sign sector.
    E01, E02, E03 : float, optional
        Pulse field amplitudes.
    phase1, phase2, phase3 : float, optional
        Carrier phases.

    Returns
    -------
    complex
        Finite-pulse third-order response.
    """
    s1, s2, s3 = pathway_signs(pathway)

    F1 = pulse_dressing_1(
        L_super,
        omega1=omega1,
        omega_L=omega_L1,
        sigma=sigma1,
        E0=E01,
        phase=phase1,
        sign=s1,
    )

    F2 = pulse_dressing_2(
        L_super,
        omega1=omega1,
        omega_L=omega_L2,
        sigma=sigma2,
        E0=E02,
        phase=phase2,
        sign=s2,
    )

    F3 = pulse_dressing_3(
        L_super,
        omega3=omega3,
        omega_L=omega_L3,
        sigma=sigma3,
        E0=E03,
        phase=phase3,
        sign=s3,
    )

    return third_order_response(
        L_super=L_super,
        V=V,
        observable=observable,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        F1=F1,
        F2=F2,
        F3=F3,
    )

def impulsive_response(
    L_super,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
):
    """
    Evaluate the impulsive third-order response using the full
    light-matter interaction superoperator.

    P^(3)(omega3, T, omega1)
        =
        <<mu|
        G(omega3) V exp(L T)
        V G(omega1) V
        |rho0>>.
    """
    L_super = np.asarray(L_super, dtype=complex)
    V = np.asarray(V, dtype=complex)
    rho0 = np.asarray(rho0, dtype=complex)

    if rho0.ndim == 2:
        state = vec(rho0)
    elif rho0.ndim == 1:
        state = rho0.copy()
    else:
        raise ValueError(
            "rho0 must be a density matrix or vectorized state"
        )

    # Pulse 1
    state = V @ state

    # Coherence-frequency interval
    state = resolvent_action(
        L_super,
        omega=omega1,
        rhs=state,
        eta=eta,
    )

    # Pulse 2
    state = V @ state

    # Waiting time
    state = expm(L_super * T) @ state

    # Pulse 3
    state = V @ state

    # Detection-frequency interval
    state = resolvent_action(
        L_super,
        omega=omega3,
        rhs=state,
        eta=eta,
    )

    bra = observable_bra(observable)

    return bra @ state
