import numpy as np

import pulsus
from dimer_model import build_dimer

def test_impulsive_signal_matches_low_level_response():
    model = build_dimer()

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        collapse_ops=collapse_ops,
        rho0=model["rho0"],
        hbar=model["hbar"],
    )

    high_level = pulsus.impulsive_signal(
        system=system,
        omega1=1.0,
        omega3=1.15,
        T=4.0,
        eta=0.02,
    )

    low_level = pulsus.impulsive_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=system.rho0,
        omega1=1.0,
        omega3=1.15,
        T=4.0,
        eta=0.02,
    )

    assert np.allclose(
        high_level,
        low_level,
    )

def test_high_level_rwa_signals_match_low_level():
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

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.50,
    )

    omega1 = 1.0
    omega3 = 1.15
    T = 4.0
    eta = 0.02

    impulsive_high = pulsus.impulsive_rwa_signal(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway="NR",
    )

    impulsive_low = pulsus.impulsive_rwa_response(
        L_super=system.L,
        V_plus=system.V_plus,
        V_minus=system.V_minus,
        observable=system.dipole,
        rho0=system.rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway="NR",
    )

    assert np.allclose(
        impulsive_high,
        impulsive_low,
    )

    short_high = pulsus.short_pulse_rwa_signal(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway="NR",
    )

    short_low = pulsus.short_pulse_rwa_response(
        L_super=system.L,
        V_plus=system.V_plus,
        V_minus=system.V_minus,
        observable=system.dipole,
        rho0=system.rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        omega_L1=pulse.omega_L,
        omega_L2=pulse.omega_L,
        omega_L3=pulse.omega_L,
        sigma1=pulse.sigma,
        sigma2=pulse.sigma,
        sigma3=pulse.sigma,
        pathway="NR",
        E01=pulse.E0,
        E02=pulse.E0,
        E03=pulse.E0,
        phase1=pulse.phase,
        phase2=pulse.phase,
        phase3=pulse.phase,
    )

    assert np.allclose(
        short_high,
        short_low,
    )

def test_rwa_signal_requires_dipole_decomposition():
    model = build_dimer()

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        rho0=model["rho0"],
    )

    try:
        pulsus.impulsive_rwa_signal(
            system=system,
            omega1=1.0,
            omega3=1.1,
            T=4.0,
            eta=0.02,
        )

    except ValueError:
        pass

    else:
        raise AssertionError(
            "RWA signal should require dipole_plus and dipole_minus"
        )

def test_finite_pulse_signal_matches_low_level_response():
    model = build_dimer()

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    system = pulsus.SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        collapse_ops=collapse_ops,
        rho0=model["rho0"],
        hbar=model["hbar"],
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

    high_level = pulsus.finite_pulse_signal(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=1.00,
        omega3=1.15,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    low_level = pulsus.finite_pulse_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=system.rho0,
        omega1=1.00,
        omega3=1.15,
        T=4.0,
        eta=0.02,
        omega_L1=pulse1.omega_L,
        omega_L2=pulse2.omega_L,
        omega_L3=pulse3.omega_L,
        sigma1=pulse1.sigma,
        sigma2=pulse2.sigma,
        sigma3=pulse3.sigma,
        pathway="NR",
        E01=pulse1.E0,
        E02=pulse2.E0,
        E03=pulse3.E0,
        phase1=pulse1.phase,
        phase2=pulse2.phase,
        phase3=pulse3.phase,
    )

    assert np.allclose(
        high_level,
        low_level,
    )

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
