import numpy as np

import pulsus

def test_arbitrary_finite_pulse_signature_matches_point_grid():
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

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.50,
    )

    omega1 = np.array([
        0.95,
        1.05,
        1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    signature = (
        1,
        1,
        1,
    )

    optimized = pulsus.finite_pulse_spectrum(
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

    reference = pulsus.response_grid(
        response_function=pulsus.finite_pulse_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        T=4.0,
        eta=0.02,
        signature=signature,
    )

    assert np.allclose(
        optimized,
        reference,
    )

def test_arbitrary_rwa_signature_matches_point_grid():
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

    # First and third field signs are negative, so use
    # negative signed frequency grids for this test.
    omega1 = np.array([
        -1.15,
        -1.05,
        -0.95,
    ])

    omega3 = np.array([
        -1.20,
        -1.00,
    ])

    signature = (
        -1,
        1,
        -1,
    )

    impulsive_fast = pulsus.impulsive_rwa_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        signature=signature,
    )

    impulsive_reference = pulsus.response_grid(
        response_function=pulsus.impulsive_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=4.0,
        eta=0.02,
        signature=signature,
    )

    assert np.allclose(
        impulsive_fast,
        impulsive_reference,
    )

    short_fast = pulsus.short_pulse_rwa_spectrum(
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

    short_reference = pulsus.response_grid(
        response_function=pulsus.short_pulse_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        T=4.0,
        eta=0.02,
        signature=signature,
    )

    assert np.allclose(
        short_fast,
        short_reference,
    )

def test_optimized_impulsive_spectrum_matches_negative_frequency_grid():
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

    omega1 = np.array([
        -0.95,
        -1.05,
        -1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    optimized = pulsus.impulsive_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
    )

    reference = pulsus.response_grid(
        response_function=pulsus.impulsive_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=4.0,
        eta=0.02,
    )

    assert optimized.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        optimized,
        reference,
    )

def test_optimized_rwa_spectra_match_response_grid_rephasing():
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

    omega1 = np.array([
        -0.95,
        -1.05,
        -1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    # --------------------------------------------------------
    # Impulsive RWA
    # --------------------------------------------------------

    impulsive_fast = pulsus.impulsive_rwa_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    impulsive_reference = pulsus.response_grid(
        response_function=pulsus.impulsive_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    assert np.allclose(
        impulsive_fast,
        impulsive_reference,
    )

    # --------------------------------------------------------
    # Short-pulse RWA
    # --------------------------------------------------------

    short_fast = pulsus.short_pulse_rwa_spectrum(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    short_reference = pulsus.response_grid(
        response_function=pulsus.short_pulse_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    assert np.allclose(
        short_fast,
        short_reference,
    )

def test_finite_pulse_spectrum_matches_response_grid_rephasing():
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

    pulse = pulsus.GaussianPulse(
        omega_L=1.10,
        sigma=0.50,
    )

    omega1 = np.array([
        -0.95,
        -1.05,
        -1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    optimized = pulsus.finite_pulse_spectrum(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    reference = pulsus.response_grid(
        response_function=pulsus.finite_pulse_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        T=4.0,
        eta=0.02,
        pathway="R",
    )

    assert optimized.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        optimized,
        reference,
    )

def test_impulsive_spectrum_matches_response_grid():
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

    omega1 = np.array([
        0.95,
        1.05,
        1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    high_level = pulsus.impulsive_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
    )

    reference = pulsus.response_grid(
        response_function=pulsus.impulsive_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=4.0,
        eta=0.02,
    )

    assert high_level.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        high_level,
        reference,
    )

def test_high_level_rwa_spectra_match_response_grid():
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

    omega1 = np.array([
        0.95,
        1.05,
        1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    # --------------------------------------------------------
    # Short-pulse RWA
    # --------------------------------------------------------

    short_high = pulsus.short_pulse_rwa_spectrum(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    short_reference = pulsus.response_grid(
        response_function=pulsus.short_pulse_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    assert short_high.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        short_high,
        short_reference,
    )

    # --------------------------------------------------------
    # Impulsive RWA
    # --------------------------------------------------------

    impulsive_high = pulsus.impulsive_rwa_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    impulsive_reference = pulsus.response_grid(
        response_function=pulsus.impulsive_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    assert impulsive_high.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        impulsive_high,
        impulsive_reference,
    )

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
def test_finite_pulse_spectrum_matches_response_grid():
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

    omega1 = np.array([
        0.95,
        1.05,
        1.15,
    ])

    omega3 = np.array([
        1.00,
        1.20,
    ])

    high_level = pulsus.finite_pulse_spectrum(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=omega1,
        omega3=omega3,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    reference = pulsus.response_grid(
        response_function=pulsus.finite_pulse_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    assert high_level.shape == (
        omega3.size,
        omega1.size,
    )

    assert np.allclose(
        high_level,
        reference,
    )
