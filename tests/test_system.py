import numpy as np

import pulsus
from dimer_model import build_dimer


def test_gaussian_pulse():
    pulse = pulsus.GaussianPulse(
        omega_L=1.1,
        sigma=1.25,
        E0=0.8,
        phase=0.2,
    )

    assert pulse.omega_L == 1.1
    assert pulse.sigma == 1.25
    assert pulse.E0 == 0.8
    assert pulse.phase == 0.2

def test_spectroscopy_system_matches_dimer_reference():
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

    assert system.d == model["d"]

    assert np.allclose(
        system.L,
        model["L_super"],
    )

    assert np.allclose(
        system.V,
        model["V"],
    )

    rho_ss = system.stationary_state()

    assert np.allclose(
        np.trace(rho_ss),
        1.0,
    )

    assert np.linalg.norm(
        system.L @ pulsus.vec(rho_ss)
    ) < 1e-10
