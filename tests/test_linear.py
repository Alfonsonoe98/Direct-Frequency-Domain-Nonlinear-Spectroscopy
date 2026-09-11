import numpy as np
import pytest

import pulsus
from dimer_model import build_dimer


def build_linear_system():
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
        hbar=model["hbar"],
        dipole_plus=mu_plus,
        dipole_minus=mu_minus,
    )

    return model, system


def test_impulsive_linear_response_matches_explicit_matrix():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    omega = 1.10
    eta = 0.02

    result = pulsus.impulsive_linear_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=rho_ss,
        omega=omega,
        eta=eta,
    )

    G = pulsus.resolvent_matrix(
        L_super=system.L,
        omega=omega,
        eta=eta,
    )

    expected = (
        pulsus.observable_bra(system.dipole)
        @ G
        @ system.V
        @ pulsus.vec(rho_ss)
    )

    assert np.allclose(
        result,
        expected,
    )


def test_finite_linear_stationary_factorization():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=1.25,
        E0=0.8,
        phase=0.2,
    )

    omega = 1.20
    eta = 0.02
    sign = 1

    finite = pulsus.finite_pulse_linear_signal(
        system=system,
        pulse=pulse,
        omega=omega,
        eta=eta,
        sign=sign,
        rho0=rho_ss,
    )

    molecular = pulsus.impulsive_linear_signal(
        system=system,
        omega=omega,
        eta=eta,
        rho0=rho_ss,
    )

    field = pulsus.gaussian_spectrum(
        omega=omega,
        omega_L=pulse.omega_L,
        sigma=pulse.sigma,
        E0=pulse.E0,
        phase=pulse.phase,
        sign=sign,
    )

    assert np.allclose(
        finite,
        field * molecular,
        rtol=1e-11,
        atol=1e-12,
    )


def test_impulsive_linear_rwa_selects_positive_interaction():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    omega = 1.10
    eta = 0.02

    result = pulsus.impulsive_linear_rwa_signal(
        system=system,
        omega=omega,
        eta=eta,
        sign=1,
        rho0=rho_ss,
    )

    G = pulsus.resolvent_matrix(
        L_super=system.L,
        omega=omega,
        eta=eta,
    )

    expected = (
        pulsus.observable_bra(system.dipole)
        @ G
        @ system.V_plus
        @ pulsus.vec(rho_ss)
    )

    assert np.allclose(
        result,
        expected,
    )


def test_impulsive_linear_rwa_selects_negative_interaction():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    omega = -1.10
    eta = 0.02

    result = pulsus.impulsive_linear_rwa_signal(
        system=system,
        omega=omega,
        eta=eta,
        sign=-1,
        rho0=rho_ss,
    )

    G = pulsus.resolvent_matrix(
        L_super=system.L,
        omega=omega,
        eta=eta,
    )

    expected = (
        pulsus.observable_bra(system.dipole)
        @ G
        @ system.V_minus
        @ pulsus.vec(rho_ss)
    )

    assert np.allclose(
        result,
        expected,
    )


def test_short_pulse_linear_rwa_factorization():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=1.25,
        E0=0.8,
        phase=0.2,
    )

    omega = 1.20
    eta = 0.02
    sign = 1

    short = pulsus.short_pulse_linear_rwa_signal(
        system=system,
        pulse=pulse,
        omega=omega,
        eta=eta,
        sign=sign,
        rho0=rho_ss,
    )

    molecular = pulsus.impulsive_linear_rwa_signal(
        system=system,
        omega=omega,
        eta=eta,
        sign=sign,
        rho0=rho_ss,
    )

    field = pulsus.gaussian_spectrum(
        omega=omega,
        omega_L=pulse.omega_L,
        sigma=pulse.sigma,
        E0=pulse.E0,
        phase=pulse.phase,
        sign=sign,
    )

    assert np.allclose(
        short,
        field * molecular,
    )


def test_linear_spectra_match_point_calls():
    model, system = build_linear_system()

    rho_ss = system.stationary_state()

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=1.25,
    )

    omega = np.array([
        0.90,
        1.00,
        1.10,
        1.20,
        1.30,
    ])

    impulsive = pulsus.impulsive_linear_spectrum(
        system=system,
        omega=omega,
        eta=0.02,
        rho0=rho_ss,
    )

    finite = pulsus.finite_pulse_linear_spectrum(
        system=system,
        pulse=pulse,
        omega=omega,
        eta=0.02,
        sign=1,
        rho0=rho_ss,
    )

    impulsive_rwa = pulsus.impulsive_linear_rwa_spectrum(
        system=system,
        omega=omega,
        eta=0.02,
        sign=1,
        rho0=rho_ss,
    )

    short_rwa = pulsus.short_pulse_linear_rwa_spectrum(
        system=system,
        pulse=pulse,
        omega=omega,
        eta=0.02,
        sign=1,
        rho0=rho_ss,
    )

    for i, w in enumerate(omega):

        assert np.allclose(
            impulsive[i],
            pulsus.impulsive_linear_signal(
                system=system,
                omega=w,
                eta=0.02,
                rho0=rho_ss,
            ),
        )

        assert np.allclose(
            finite[i],
            pulsus.finite_pulse_linear_signal(
                system=system,
                pulse=pulse,
                omega=w,
                eta=0.02,
                sign=1,
                rho0=rho_ss,
            ),
        )

        assert np.allclose(
            impulsive_rwa[i],
            pulsus.impulsive_linear_rwa_signal(
                system=system,
                omega=w,
                eta=0.02,
                sign=1,
                rho0=rho_ss,
            ),
        )

        assert np.allclose(
            short_rwa[i],
            pulsus.short_pulse_linear_rwa_signal(
                system=system,
                pulse=pulse,
                omega=w,
                eta=0.02,
                sign=1,
                rho0=rho_ss,
            ),
        )


def test_linear_rwa_requires_dipole_decomposition():
    model = build_dimer()

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        rho0=model["rho0"],
    )

    with pytest.raises(ValueError):
        pulsus.impulsive_linear_rwa_signal(
            system=system,
            omega=1.1,
            eta=0.02,
            sign=1,
        )
