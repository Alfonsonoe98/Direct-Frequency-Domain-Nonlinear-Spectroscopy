"""A--D models with explicit NR/R axes and unit Gaussian prefactors (phase=0)."""

import matplotlib.pyplot as plt
import numpy as np
import pulsus

from common import build_dimer


def main():
    system, rho_ss = build_dimer()
    sigma = 1.25
    # E0 fixes the signed-envelope prefactor to one, not an oscillatory-field integral.
    pulse = pulsus.GaussianPulse(
        omega_L=1.10, sigma=sigma, E0=2 / (np.sqrt(2 * np.pi) * sigma)
    )
    omega_positive = np.linspace(0.7, 1.5, 41)
    models = {
        "A: finite, full": pulsus.finite_pulse_spectrum,
        "B: finite, RWA": pulsus.finite_pulse_rwa_spectrum,
        "C: short, RWA": pulsus.short_pulse_rwa_spectrum,
        "D: impulsive, RWA": pulsus.impulsive_rwa_spectrum,
    }
    fig, axes = plt.subplots(2, 4, figsize=(13, 6), constrained_layout=True)
    for row, pathway in enumerate(["NR", "R"]):
        omega1 = omega_positive if pathway == "NR" else -omega_positive[::-1]
        shared = dict(system=system, omega1=omega1, omega3=omega_positive,
                      T=20.0, eta=0.02, pathway=pathway, rho0=rho_ss)
        results = []
        for col, (name, function) in enumerate(models.items()):
            if name.startswith("D:"):
                spectrum = function(**shared)
            else:
                spectrum = function(pulse1=pulse, pulse2=pulse, pulse3=pulse, **shared)
            assert spectrum.shape == (41, 41) and np.isfinite(spectrum).all()
            results.append(spectrum)
            # Only the displayed R coordinate is sign-flipped.
            display_axis = omega1 if pathway == "NR" else -omega1
            pulsus.plot_spectrum(
                display_axis, omega_positive, spectrum, component="imag",
                normalize=True, ax=axes[row, col], colorbar=False,
                title=f"{pathway}: {name}",
                xlabel=r"$\omega_1$" if pathway == "NR" else r"$-\omega_1$",
                ylabel=r"$\omega_3$",
            )
        changes = [np.linalg.norm(x - y) / np.linalg.norm(x)
                   for x, y in zip(results[:-1], results[1:])]
        print(pathway, "raw complex relative changes A->B, B->C, C->D:", changes)
    fig.suptitle("Model units; each panel normalized separately (imaginary component)")
    print("NR and R grids:", (41, 41), "; four models per sector")
    plt.show()


if __name__ == "__main__":
    main()
