import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pulsus


def test_plot_spectrum():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.array([
        [1.0 + 0.2j, -0.5 + 0.1j, 0.2 - 0.1j],
        [-0.8 + 0.3j, 0.4 - 0.2j, 0.1 + 0.5j],
    ])

    fig, ax, contour = pulsus.plot_spectrum(
        omega1=omega1,
        omega3=omega3,
        spectrum=spectrum,
        component="real",
    )

    assert fig is not None
    assert ax is not None
    assert contour is not None

    assert ax.get_xlabel() == r"$\omega_1$"
    assert ax.get_ylabel() == r"$\omega_3$"

    plt.close(fig)


def test_plot_spectrum_accepts_components():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.ones(
        (2, 3),
        dtype=complex,
    )

    for component in (
        "real",
        "imag",
        "abs",
    ):
        fig, ax, contour = pulsus.plot_spectrum(
            omega1=omega1,
            omega3=omega3,
            spectrum=spectrum,
            component=component,
            normalize=True,
        )

        assert contour is not None

        plt.close(fig)


def test_plot_spectrum_rejects_wrong_shape():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    wrong_shape = np.zeros(
        (3, 2),
        dtype=complex,
    )

    with pytest.raises(ValueError):
        pulsus.plot_spectrum(
            omega1=omega1,
            omega3=omega3,
            spectrum=wrong_shape,
        )
