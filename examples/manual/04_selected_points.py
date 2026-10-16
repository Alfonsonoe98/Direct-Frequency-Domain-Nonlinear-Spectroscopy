"""Paired frequency points checked against a rectangular Cartesian grid."""

import matplotlib.pyplot as plt
import numpy as np
import pulsus

from common import build_dimer


def main():
    system, rho_ss = build_dimer()
    pulse = pulsus.GaussianPulse(omega_L=1.10, sigma=1.25)
    omega1 = np.linspace(0.7, 1.5, 17)
    omega3 = np.linspace(0.7, 1.5, 13)
    shared = dict(system=system, pulse1=pulse, pulse2=pulse, pulse3=pulse,
                  T=20.0, eta=0.02, pathway="NR", rho0=rho_ss)
    grid = pulsus.finite_pulse_spectrum(omega1=omega1, omega3=omega3, **shared)
    indices = np.array([[2, 3], [6, 8], [8, 11], [10, 14]])  # (row j, column i)
    points1 = omega1[indices[:, 1]]
    points3 = omega3[indices[:, 0]]
    values = pulsus.finite_pulse_points(omega1=points1, omega3=points3, **shared)
    expected = grid[indices[:, 0], indices[:, 1]]
    np.testing.assert_allclose(values, expected, rtol=1e-10, atol=1e-10)
    assert grid.shape == (13, 17) and values.shape == (4,)
    assert np.isfinite(grid).all() and np.isfinite(values).all()
    print("Cartesian grid shape:", grid.shape)
    print("Paired points shape:", values.shape)
    print("Maximum points/grid difference:", np.max(np.abs(values - expected)))
    fig, ax, _ = pulsus.plot_spectrum(
        omega1, omega3, grid, component="imag", normalize=True,
        title="Four selected coordinates on a rectangular NR grid",
        xlabel=r"$\omega_1$ (inverse model-time units)",
        ylabel=r"$\omega_3$ (inverse model-time units)",
    )
    ax.scatter(points1, points3, facecolors="none", edgecolors="black", s=70)
    for k, (x, y) in enumerate(zip(points1, points3), start=1):
        ax.annotate(str(k), (x, y), xytext=(5, 5), textcoords="offset points")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
