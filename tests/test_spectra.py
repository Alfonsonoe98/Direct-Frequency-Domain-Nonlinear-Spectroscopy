import numpy as np

import pulsus


def test_response_grid():
    omega1 = np.array([0.9, 1.0, 1.1])
    omega3 = np.array([1.0, 1.2])

    def mock_response(omega1, omega3, scale):
        return scale * (omega1 + 1j * omega3)

    spectrum = pulsus.response_grid(
        response_function=mock_response,
        omega1=omega1,
        omega3=omega3,
        scale=2.0,
    )

    assert spectrum.shape == (2, 3)

    for j, w3 in enumerate(omega3):
        for i, w1 in enumerate(omega1):
            reference = 2.0 * (w1 + 1j * w3)

            assert np.allclose(
                spectrum[j, i],
                reference,
            )

from dimer_model import build_dimer


def test_response_grid_matches_direct_impulsive_calls():
    model = build_dimer()

    omega1 = np.array([0.95, 1.05, 1.15])
    omega3 = np.array([1.00, 1.20])

    spectrum = pulsus.response_grid(
        response_function=pulsus.impulsive_response,
        omega1=omega1,
        omega3=omega3,
        L_super=model["L_super"],
        V=model["V"],
        observable=model["mu"],
        rho0=model["rho0"],
        T=4.0,
        eta=0.02,
    )

    for j, w3 in enumerate(omega3):
        for i, w1 in enumerate(omega1):

            direct = pulsus.impulsive_response(
                L_super=model["L_super"],
                V=model["V"],
                observable=model["mu"],
                rho0=model["rho0"],
                omega1=w1,
                omega3=w3,
                T=4.0,
                eta=0.02,
            )

            assert np.allclose(
                spectrum[j, i],
                direct,
            )
