import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import expm

import pulsus
from dimer_model import build_dimer

def _gaussian_dressing_integral(
    L_super,
    detuning,
    sigma,
    E0,
    phase,
    sign,
):
    """
    Direct numerical evaluation of the defining Gaussian
    time-domain dressing integral.
    """
    L_super = np.asarray(L_super, dtype=complex)

    n = L_super.shape[0]
    I = np.eye(n, dtype=complex)

    A = L_super + 1j * detuning * I

    def integrand(t):
        envelope = (
            0.5
            * E0
            * np.exp(-0.5 * (t / sigma) ** 2)
            * np.exp(1j * sign * phase)
        )

        return envelope * expm(A * t)

    # ±8 sigma captures essentially the entire Gaussian.
    result, error = quad_vec(
        integrand,
        -8.0 * sigma,
        8.0 * sigma,
        epsabs=1e-12,
        epsrel=1e-12,
    )

    return result

def test_gaussian_prefactor():
    sigma = 0.8
    E0 = 1.3
    phase = 0.4
    sign = -1

    result = pulsus.gaussian_prefactor(
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )

    reference = (
        0.5
        * E0
        * np.sqrt(2.0 * np.pi)
        * sigma
        * np.exp(1j * sign * phase)
    )

    assert np.allclose(result, reference)


def test_first_pulse_stationary_state_identity():
    model = build_dimer()

    L = model["L_super"]

    rho_ss = pulsus.stationary_state(
        L,
        model["d"],
    )

    rho_ss_vec = pulsus.vec(rho_ss)

    omega1 = 1.0
    omega_L = 1.1
    sigma = 0.75
    E0 = 0.8
    phase = 0.3
    sign = 1

    F1 = pulsus.pulse_dressing_1(
        L,
        omega1=omega1,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )

    C = (
        0.5
        * E0
        * np.sqrt(2.0 * np.pi)
        * sigma
    )

    E_omega = (
        C
        * np.exp(1j * sign * phase)
        * np.exp(
            -0.5
            * sigma**2
            * (omega1 - sign * omega_L)**2
        )
    )

    lhs = F1 @ rho_ss_vec
    rhs = E_omega * rho_ss_vec

    assert np.allclose(lhs, rhs)


def test_pulse_dressing_matches_defining_integral():
    model = build_dimer()

    L = model["L_super"]

    sigma = 0.45
    E0 = 0.8
    phase = 0.23

    omega1 = 1.03
    omega3 = 1.17
    omega_L = 1.10

    # --------------------------------------------------------
    # Pulse 1
    # --------------------------------------------------------

    sign1 = 1

    F1_closed = pulsus.pulse_dressing_1(
        L,
        omega1=omega1,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign1,
    )

    detuning1 = omega1 - sign1 * omega_L

    F1_integral = _gaussian_dressing_integral(
        L,
        detuning=detuning1,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign1,
    )

    # --------------------------------------------------------
    # Pulse 2
    # --------------------------------------------------------

    sign2 = 1

    F2_closed = pulsus.pulse_dressing_2(
        L,
        omega1=omega1,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign2,
    )

    detuning2 = omega1 + sign2 * omega_L

    F2_integral = _gaussian_dressing_integral(
        L,
        detuning=detuning2,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign2,
    )

    # --------------------------------------------------------
    # Pulse 3
    # --------------------------------------------------------

    sign3 = -1

    F3_closed = pulsus.pulse_dressing_3(
        L,
        omega3=omega3,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign3,
    )

    detuning3 = omega3 - sign3 * omega_L

    F3_integral = _gaussian_dressing_integral(
        L,
        detuning=detuning3,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign3,
    )

    # --------------------------------------------------------
    # Compare
    # --------------------------------------------------------

    assert np.allclose(
        F1_closed,
        F1_integral,
        atol=1e-11,
        rtol=1e-11,
    )

    assert np.allclose(
        F2_closed,
        F2_integral,
        atol=1e-11,
        rtol=1e-11,
    )

    assert np.allclose(
        F3_closed,
        F3_integral,
        atol=1e-11,
        rtol=1e-11,
    )
