import numpy as np
from scipy.linalg import expm

import pulsus
from dimer_model import build_dimer


def test_finite_pulse_rwa_matches_manual_construction():
    model = build_dimer()

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    mu_plus = (
        model["sigma_a_plus"]
        + model["sigma_b_plus"]
    )

    mu_minus = (
        model["sigma_a_minus"]
        + model["sigma_b_minus"]
    )

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        collapse_ops=collapse_ops,
        rho0=model["rho0"],
        hbar=model["hbar"],
        dipole_plus=mu_plus,
        dipole_minus=mu_minus,
    )

    pulse1 = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.45,
    )

    pulse2 = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.50,
    )

    pulse3 = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.55,
    )

    omega1 = 1.00
    omega3 = 1.15
    T = 4.0
    eta = 0.02

    signature = (
        +1,
        -1,
        +1,
    )

    s1, s2, s3 = signature

    V1 = (
        system.V_plus
        if s1 == 1
        else system.V_minus
    )

    V2 = (
        system.V_plus
        if s2 == 1
        else system.V_minus
    )

    V3 = (
        system.V_plus
        if s3 == 1
        else system.V_minus
    )

    F1 = pulsus.pulse_dressing_1(
        system.L,
        omega1,
        pulse1.omega_L,
        pulse1.sigma,
        E0=pulse1.E0,
        phase=pulse1.phase,
        sign=s1,
    )

    F2 = pulsus.pulse_dressing_2(
        system.L,
        omega1,
        pulse2.omega_L,
        pulse2.sigma,
        E0=pulse2.E0,
        phase=pulse2.phase,
        sign=s2,
    )

    F3 = pulsus.pulse_dressing_3(
        system.L,
        omega3,
        pulse3.omega_L,
        pulse3.sigma,
        E0=pulse3.E0,
        phase=pulse3.phase,
        sign=s3,
    )

    state = pulsus.vec(
        system.rho0
    )

    state = F1 @ state
    state = V1 @ state

    state = pulsus.resolvent_action(
        L_super=system.L,
        omega=omega1,
        rhs=state,
        eta=eta,
    )

    state = V2 @ state

    state = F2 @ state
    state = F3 @ state

    state = (
        expm(system.L * T)
        @ state
    )

    state = V3 @ state

    state = pulsus.resolvent_action(
        L_super=system.L,
        omega=omega3,
        rhs=state,
        eta=eta,
    )

    manual = (
        pulsus.observable_bra(
            system.dipole
        )
        @ state
    )

    automatic = pulsus.finite_pulse_rwa_signal(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        signature=signature,
    )

    assert np.allclose(
        automatic,
        manual,
        rtol=1e-12,
        atol=1e-12,
    )
