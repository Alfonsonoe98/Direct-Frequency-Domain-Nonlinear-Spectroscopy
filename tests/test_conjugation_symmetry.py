import numpy as np

import pulsus

from dimer_model import (
    build_dimer,
    stationary_state,
)


def test_finite_pulse_sign_reversed_conjugation_symmetry():
    model = build_dimer()

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    rho_ss = stationary_state(
        model["L_super"],
        model["d"],
    )

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        collapse_ops=collapse_ops,
        rho0=rho_ss,
        hbar=model["hbar"],
    )

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.50,
    )

    omega1 = np.array([
        0.95,
        1.10,
        1.25,
    ])

    omega3 = np.array([
        0.95,
        1.10,
        1.25,
    ])

    signatures = [
        (+1, +1, +1),
        (+1, -1, +1),
        (-1, +1, +1),
        (-1, -1, +1),
    ]

    for signature in signatures:
        reversed_signature = tuple(
            -sign
            for sign in signature
        )

        spectrum = pulsus.finite_pulse_spectrum(
            system=system,
            pulse1=pulse,
            pulse2=pulse,
            pulse3=pulse,
            omega1=omega1,
            omega3=omega3,
            T=4.0,
            eta=0.02,
            signature=signature,
        )

        reversed_spectrum = pulsus.finite_pulse_spectrum(
            system=system,
            pulse1=pulse,
            pulse2=pulse,
            pulse3=pulse,
            omega1=-omega1,
            omega3=-omega3,
            T=4.0,
            eta=0.02,
            signature=reversed_signature,
        )

        np.testing.assert_allclose(
            reversed_spectrum,
            spectrum.conj(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
