import numpy as np
from scipy.linalg import expm

from dimer_model import build_dimer

from pulsus.system import SpectroscopySystem
from pulsus.liouville import vec
from pulsus.response import observable_bra
from pulsus.time_domain import (
    impulsive_linear_time_signal,
    impulsive_linear_rwa_time_signal,
    impulsive_third_order_time_signal,
    impulsive_third_order_rwa_time_signal,
    impulsive_linear_polarization,
    impulsive_third_order_polarization,
    impulsive_pathway_polarization,
)

def test_impulsive_pathway_polarization_piecewise_history():
    system = build_system()

    pulse_times = (
        1.0,
        3.0,
        7.0,
    )

    times = np.array([
        0.0,
        1.0,
        2.0,
        3.0,
        5.0,
        7.0,
        8.0,
    ])

    signature = (
        +1,
        -1,
        +1,
    )

    polarization, components = (
        impulsive_pathway_polarization(
            system=system,
            times=times,
            pulse_times=pulse_times,
            signature=signature,
            return_components=True,
        )
    )

    P1 = components["first_order"]
    P2 = components["second_order"]
    P3 = components["third_order"]

    assert np.allclose(
        polarization[0],
        0.0,
    )

    # First-order stage
    expected1 = (
        observable_bra(system.dipole)
        @ expm(system.L * 1.0)
        @ system.V_plus
        @ vec(system.rho0)
    )

    assert np.allclose(
        P1[2],
        expected1,
    )

    # Second-order stage
    expected2 = (
        observable_bra(system.dipole)
        @ expm(system.L * 2.0)
        @ system.V_minus
        @ expm(system.L * 2.0)
        @ system.V_plus
        @ vec(system.rho0)
    )

    assert np.allclose(
        P2[4],
        expected2,
    )

    # Third-order stage
    expected3 = (
        observable_bra(system.dipole)
        @ expm(system.L * 1.0)
        @ system.V_plus
        @ expm(system.L * 4.0)
        @ system.V_minus
        @ expm(system.L * 2.0)
        @ system.V_plus
        @ vec(system.rho0)
    )

    assert np.allclose(
        P3[6],
        expected3,
    )

    assert np.allclose(
        polarization,
        P1 + P2 + P3,
    )

def build_system():
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

    system = SpectroscopySystem(
        H=model["H_S"],
        dipole=model["mu"],
        collapse_ops=collapse_ops,
        rho0=model["rho0"],
        hbar=model["hbar"],
        dipole_plus=mu_plus,
        dipole_minus=mu_minus,
    )

    return system


def test_impulsive_linear_time_zero():
    system = build_system()

    signal = impulsive_linear_time_signal(
        system=system,
        times=np.array([0.0]),
    )

    expected = (
        observable_bra(system.dipole)
        @ system.V
        @ vec(system.rho0)
    )

    assert np.allclose(
        signal[0],
        expected,
    )


def test_impulsive_linear_rwa_time_zero():
    system = build_system()

    signal = impulsive_linear_rwa_time_signal(
        system=system,
        times=np.array([0.0]),
        sign=1,
    )

    expected = (
        observable_bra(system.dipole)
        @ system.V_plus
        @ vec(system.rho0)
    )

    assert np.allclose(
        signal[0],
        expected,
    )


def test_impulsive_third_order_time_matches_manual_point():
    system = build_system()

    t1 = np.array([
        0.3,
        0.8,
    ])

    t3 = np.array([
        0.2,
        0.6,
    ])

    T = 4.0

    signal = impulsive_third_order_time_signal(
        system=system,
        t1=t1,
        t3=t3,
        T=T,
    )

    expected = (
        observable_bra(system.dipole)
        @ expm(system.L * t3[1])
        @ system.V
        @ expm(system.L * T)
        @ system.V
        @ expm(system.L * t1[0])
        @ system.V
        @ vec(system.rho0)
    )

    assert np.allclose(
        signal[1, 0],
        expected,
    )


def test_impulsive_third_order_rwa_time_matches_manual_point():
    system = build_system()

    t1 = np.array([
        0.3,
        0.8,
    ])

    t3 = np.array([
        0.2,
        0.6,
    ])

    T = 4.0

    signature = (
        -1,
        +1,
        +1,
    )

    signal = impulsive_third_order_rwa_time_signal(
        system=system,
        t1=t1,
        t3=t3,
        T=T,
        signature=signature,
    )

    expected = (
        observable_bra(system.dipole)
        @ expm(system.L * t3[1])
        @ system.V_plus
        @ expm(system.L * T)
        @ system.V_plus
        @ expm(system.L * t1[0])
        @ system.V_minus
        @ vec(system.rho0)
    )

    assert np.allclose(
        signal[1, 0],
        expected,
    )


def test_impulsive_linear_polarization_matches_shifted_response():
    system = build_system()

    times = np.array([
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
    ])

    pulse_time = 2.0

    polarization = impulsive_linear_polarization(
        system=system,
        times=times,
        pulse_time=pulse_time,
    )

    assert np.allclose(
        polarization[:2],
        0.0,
    )

    reference = impulsive_linear_time_signal(
        system=system,
        times=np.array([
            0.0,
            1.0,
            2.0,
        ]),
    )

    assert np.allclose(
        polarization[2:],
        reference,
    )


def test_impulsive_linear_rwa_polarization_matches_shifted_response():
    system = build_system()

    times = np.array([
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
    ])

    pulse_time = 2.0

    polarization = impulsive_linear_polarization(
        system=system,
        times=times,
        pulse_time=pulse_time,
        sign=+1,
    )

    reference = impulsive_linear_rwa_time_signal(
        system=system,
        times=np.array([
            0.0,
            1.0,
            2.0,
        ]),
        sign=+1,
    )

    assert np.allclose(
        polarization[:2],
        0.0,
    )

    assert np.allclose(
        polarization[2:],
        reference,
    )


def test_impulsive_third_order_polarization_matches_time_response():
    system = build_system()

    pulse_times = (
        1.0,
        3.0,
        7.0,
    )

    times = np.array([
        0.0,
        2.0,
        6.0,
        7.0,
        8.0,
        9.0,
    ])

    polarization = impulsive_third_order_polarization(
        system=system,
        times=times,
        pulse_times=pulse_times,
    )

    assert np.allclose(
        polarization[:3],
        0.0,
    )

    reference = impulsive_third_order_time_signal(
        system=system,
        t1=np.array([2.0]),
        t3=np.array([
            0.0,
            1.0,
            2.0,
        ]),
        T=4.0,
    )

    assert np.allclose(
        polarization[3:],
        reference[:, 0],
    )


def test_impulsive_third_order_rwa_polarization_matches_time_response():
    system = build_system()

    pulse_times = (
        1.0,
        3.0,
        7.0,
    )

    times = np.array([
        0.0,
        2.0,
        6.0,
        7.0,
        8.0,
        9.0,
    ])

    signature = (
        +1,
        -1,
        +1,
    )

    polarization = impulsive_third_order_polarization(
        system=system,
        times=times,
        pulse_times=pulse_times,
        signature=signature,
    )

    reference = impulsive_third_order_rwa_time_signal(
        system=system,
        t1=np.array([2.0]),
        t3=np.array([
            0.0,
            1.0,
            2.0,
        ]),
        T=4.0,
        signature=signature,
    )

    assert np.allclose(
        polarization[:3],
        0.0,
    )

    assert np.allclose(
        polarization[3:],
        reference[:, 0],
    )
