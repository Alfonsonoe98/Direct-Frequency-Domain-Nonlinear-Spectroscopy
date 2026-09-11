import numpy as np
from scipy.linalg import expm

import pulsus
from dimer_model import build_dimer

def test_resolve_signature_rejects_invalid_signatures():
    import pytest

    with pytest.raises(ValueError):
        pulsus.resolve_signature(
            signature=(1, -1),
        )

    with pytest.raises(ValueError):
        pulsus.resolve_signature(
            signature=(1, 0, 1),
        )

def test_resolve_signature():
    assert pulsus.resolve_signature(
        pathway="NR"
    ) == (
        1,
        -1,
        1,
    )

    assert pulsus.resolve_signature(
        pathway="R"
    ) == (
        -1,
        1,
        1,
    )

    assert pulsus.resolve_signature(
        pathway="NR",
        signature=(1, 1, 1),
    ) == (
        1,
        1,
        1,
    )

def test_impulsive_rwa_pathway_selection():
    model = build_dimer()

    mu_plus = (
        model["sigma_a_plus"]
        + model["sigma_b_plus"]
    )

    mu_minus = (
        model["sigma_a_minus"]
        + model["sigma_b_minus"]
    )

    V_plus = pulsus.interaction_superoperator(mu_plus)
    V_minus = pulsus.interaction_superoperator(mu_minus)

    NR = pulsus.impulsive_rwa_response(
        L_super=model["L_super"],
        V_plus=V_plus,
        V_minus=V_minus,
        observable=model["mu"],
        rho0=model["rho0"],
        omega1=1.0,
        omega3=1.15,
        T=4.0,
        eta=0.02,
        pathway="NR",
    )

    assert np.isfinite(NR.real)
    assert np.isfinite(NR.imag)

def test_short_pulse_rwa_factorization():
    model = build_dimer()

    mu_plus = (
        model["sigma_a_plus"]
        + model["sigma_b_plus"]
    )

    mu_minus = (
        model["sigma_a_minus"]
        + model["sigma_b_minus"]
    )

    V_plus = pulsus.interaction_superoperator(mu_plus)
    V_minus = pulsus.interaction_superoperator(mu_minus)

    omega1 = 1.0
    omega3 = 1.15
    T = 4.0
    eta = 0.02

    sigma = 0.5
    omega_L = 1.1

    molecular = pulsus.impulsive_rwa_response(
        L_super=model["L_super"],
        V_plus=V_plus,
        V_minus=V_minus,
        observable=model["mu"],
        rho0=model["rho0"],
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        pathway="NR",
    )

    E1 = pulsus.gaussian_spectrum(
        omega1,
        omega_L,
        sigma,
        sign=1,
    )

    E2 = pulsus.gaussian_spectrum(
        -omega1,
        omega_L,
        sigma,
        sign=-1,
    )

    E3 = pulsus.gaussian_spectrum(
        omega3,
        omega_L,
        sigma,
        sign=1,
    )

    reference = E3 * E2 * E1 * molecular

    result = pulsus.short_pulse_rwa_response(
        L_super=model["L_super"],
        V_plus=V_plus,
        V_minus=V_minus,
        observable=model["mu"],
        rho0=model["rho0"],
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        omega_L1=omega_L,
        omega_L2=omega_L,
        omega_L3=omega_L,
        sigma1=sigma,
        sigma2=sigma,
        sigma3=sigma,
        pathway="NR",
    )

    assert np.allclose(result, reference)

def test_impulsive_response_matches_identity_dressing():
    model = build_dimer()

    L = model["L_super"]
    V = model["V"]
    mu = model["mu"]
    rho0 = model["rho0"]

    omega1 = 1.0
    omega3 = 1.15
    T = 4.0
    eta = 0.02

    impulsive = pulsus.impulsive_response(
        L_super=L,
        V=V,
        observable=mu,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
    )

    I = np.eye(L.shape[0], dtype=complex)

    identity_dressed = pulsus.third_order_response(
        L_super=L,
        V=V,
        observable=mu,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        F1=I,
        F2=I,
        F3=I,
    )

    assert np.allclose(
        impulsive,
        identity_dressed,
    )

def test_pathway_signs():
    assert pulsus.pathway_signs("NR") == (1, -1, 1)
    assert pulsus.pathway_signs("R") == (-1, 1, 1)

    # Case-insensitive convenience
    assert pulsus.pathway_signs("nr") == (1, -1, 1)
    assert pulsus.pathway_signs("r") == (-1, 1, 1)

def test_finite_pulse_response_matches_manual_construction():
    model = build_dimer()

    L = model["L_super"]
    V = model["V"]
    mu = model["mu"]
    rho0 = model["rho0"]

    omega1 = 1.00
    omega3 = 1.15

    T = 4.0
    eta = 0.02

    omega_L1 = 1.10
    omega_L2 = 1.10
    omega_L3 = 1.10

    sigma1 = 0.45
    sigma2 = 0.50
    sigma3 = 0.55

    s1, s2, s3 = pulsus.pathway_signs("NR")

    F1 = pulsus.pulse_dressing_1(
        L,
        omega1,
        omega_L1,
        sigma1,
        sign=s1,
    )

    F2 = pulsus.pulse_dressing_2(
        L,
        omega1,
        omega_L2,
        sigma2,
        sign=s2,
    )

    F3 = pulsus.pulse_dressing_3(
        L,
        omega3,
        omega_L3,
        sigma3,
        sign=s3,
    )

    manual = pulsus.third_order_response(
        L_super=L,
        V=V,
        observable=mu,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        F1=F1,
        F2=F2,
        F3=F3,
    )

    automatic = pulsus.finite_pulse_response(
        L_super=L,
        V=V,
        observable=mu,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        omega_L1=omega_L1,
        omega_L2=omega_L2,
        omega_L3=omega_L3,
        sigma1=sigma1,
        sigma2=sigma2,
        sigma3=sigma3,
        pathway="NR",
    )

    assert np.allclose(automatic, manual)

def test_observable_bra():
    model = build_dimer()

    mu = model["mu"]
    rho = model["rho0"]

    bra = pulsus.observable_bra(mu)

    result = bra @ pulsus.vec(rho)
    reference = np.trace(mu @ rho)

    assert np.allclose(result, reference)


def test_third_order_response_matches_explicit_matrix_expression():
    model = build_dimer()

    L = model["L_super"]
    V = model["V"]
    mu = model["mu"]
    rho0 = model["rho0"]

    omega1 = 1.00
    omega3 = 1.15
    T = 4.0
    eta = 0.02

    sigma1 = 0.45
    sigma2 = 0.50
    sigma3 = 0.55

    omega_L1 = 1.10
    omega_L2 = 1.10
    omega_L3 = 1.10

    s1 = 1
    s2 = -1
    s3 = 1

    F1 = pulsus.pulse_dressing_1(
        L,
        omega1,
        omega_L1,
        sigma1,
        sign=s1,
    )

    F2 = pulsus.pulse_dressing_2(
        L,
        omega1,
        omega_L2,
        sigma2,
        sign=s2,
    )

    F3 = pulsus.pulse_dressing_3(
        L,
        omega3,
        omega_L3,
        sigma3,
        sign=s3,
    )

    result = pulsus.third_order_response(
        L_super=L,
        V=V,
        observable=mu,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        F1=F1,
        F2=F2,
        F3=F3,
    )

    G1 = pulsus.resolvent_matrix(
        L,
        omega1,
        eta,
    )

    G3 = pulsus.resolvent_matrix(
        L,
        omega3,
        eta,
    )

    U_T = expm(L * T)

    rho_vec = pulsus.vec(rho0)
    mu_bra = pulsus.observable_bra(mu)

    reference = (
        mu_bra
        @ G3
        @ V
        @ U_T
        @ F3
        @ F2
        @ V
        @ G1
        @ V
        @ F1
        @ rho_vec
    )

    assert np.allclose(result, reference)
