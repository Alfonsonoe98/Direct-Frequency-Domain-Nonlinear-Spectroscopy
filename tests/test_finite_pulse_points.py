import numpy as np
import pytest

import pulsus


def make_test_system():
    """Small two-level system for paired-frequency API tests."""

    H = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=complex,
    )

    dipole = np.array(
        [
            [0.0, 1.0],
            [1.0, 0.0],
        ],
        dtype=complex,
    )

    rho0 = np.array(
        [
            [1.0, 0.0],
            [0.0, 0.0],
        ],
        dtype=complex,
    )

    system = pulsus.SpectroscopySystem(
        H=H,
        dipole=dipole,
        collapse_ops=[],
        rho0=rho0,
        hbar=1.0,
    )

    pulse = pulsus.GaussianPulse(
        omega_L=1.0,
        sigma=0.5,
        E0=1.0,
        phase=0.0,
    )

    return system, pulse


def test_finite_pulse_points_matches_cartesian_grid():
    """
    Paired-coordinate evaluation must reproduce selected points
    from the ordinary Cartesian spectrum.
    """

    system, pulse = make_test_system()

    omega1_axis = np.array(
        [0.72, 0.91, 1.17, 1.34]
    )

    omega3_axis = np.array(
        [0.76, 1.03, 1.29]
    )

    grid = pulsus.finite_pulse_spectrum(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1_axis,
        omega3=omega3_axis,
        T=0.4,
        eta=1.0e-3,
        pathway="NR",
    )

    i_indices = np.array(
        [0, 3, 1, 2, 3]
    )

    j_indices = np.array(
        [2, 0, 1, 2, 1]
    )

    omega1_points = omega1_axis[
        i_indices
    ]

    omega3_points = omega3_axis[
        j_indices
    ]

    points = pulsus.finite_pulse_points(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1_points,
        omega3=omega3_points,
        T=0.4,
        eta=1.0e-3,
        pathway="NR",
    )

    expected = grid[
        j_indices,
        i_indices,
    ]

    np.testing.assert_allclose(
        points,
        expected,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_finite_pulse_points_handles_repeated_frequencies():
    """
    Repeated omega1 and omega3 values must give the same result
    as independent Cartesian-grid evaluation.
    """

    system, pulse = make_test_system()

    omega1_points = np.array(
        [
            0.85,
            0.85,
            1.15,
            1.15,
            0.85,
        ]
    )

    omega3_points = np.array(
        [
            0.90,
            1.20,
            0.90,
            1.20,
            0.90,
        ]
    )

    points = pulsus.finite_pulse_points(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1_points,
        omega3=omega3_points,
        T=0.25,
        eta=1.0e-3,
        pathway="NR",
    )

    omega1_unique = np.unique(
        omega1_points
    )

    omega3_unique = np.unique(
        omega3_points
    )

    grid = pulsus.finite_pulse_spectrum(
        system=system,
        pulse1=pulse,
        pulse2=pulse,
        pulse3=pulse,
        omega1=omega1_unique,
        omega3=omega3_unique,
        T=0.25,
        eta=1.0e-3,
        pathway="NR",
    )

    expected = np.empty(
        omega1_points.size,
        dtype=complex,
    )

    for k, (w1, w3) in enumerate(
        zip(
            omega1_points,
            omega3_points,
        )
    ):

        i = np.where(
            omega1_unique == w1
        )[0][0]

        j = np.where(
            omega3_unique == w3
        )[0][0]

        expected[k] = grid[
            j,
            i,
        ]

    np.testing.assert_allclose(
        points,
        expected,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_finite_pulse_points_requires_matching_1d_arrays():
    """omega1 and omega3 must be paired one-dimensional arrays."""

    system, pulse = make_test_system()

    with pytest.raises(
        ValueError,
        match="same shape",
    ):

        pulsus.finite_pulse_points(
            system=system,
            pulse1=pulse,
            pulse2=pulse,
            pulse3=pulse,
            omega1=np.array(
                [0.8, 1.0, 1.2]
            ),
            omega3=np.array(
                [0.9, 1.1]
            ),
            T=0.2,
            eta=1.0e-3,
            pathway="NR",
        )

    with pytest.raises(
        ValueError,
        match="one-dimensional",
    ):

        pulsus.finite_pulse_points(
            system=system,
            pulse1=pulse,
            pulse2=pulse,
            pulse3=pulse,
            omega1=np.array(
                [
                    [0.8, 1.0],
                ]
            ),
            omega3=np.array(
                [
                    [0.9, 1.1],
                ]
            ),
            T=0.2,
            eta=1.0e-3,
            pathway="NR",
        )
