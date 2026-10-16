"""Four public linear models; curves normalized separately for shape comparison."""

import matplotlib.pyplot as plt
import numpy as np
import pulsus

from common import build_dimer


def main():
    system, rho_ss = build_dimer()
    omega = np.linspace(0.5, 1.7, 241)
    sigma = 1.25
    # This E0 makes gaussian_prefactor=1 at phase=0 (unit signed-envelope area).
    pulse = pulsus.GaussianPulse(
        omega_L=1.10, sigma=sigma, E0=2 / (np.sqrt(2 * np.pi) * sigma)
    )
    shared = dict(system=system, omega=omega, eta=0.02, rho0=rho_ss)
    spectra = {
        "Finite pulse, full": pulsus.finite_pulse_linear_spectrum(
            pulse=pulse, sign=1, **shared),
        "Impulsive, full": pulsus.impulsive_linear_spectrum(**shared),
        "Short pulse, RWA": pulsus.short_pulse_linear_rwa_spectrum(
            pulse=pulse, sign=1, **shared),
        "Impulsive, RWA": pulsus.impulsive_linear_rwa_spectrum(sign=1, **shared),
    }
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    styles = ["-", "--", "-.", ":"]
    markers = ["o", "s", "^", "D"]
    for k, ((name, spectrum), style, marker) in enumerate(
            zip(spectra.items(), styles, markers)):
        assert spectrum.shape == omega.shape and np.isfinite(spectrum).all()
        scale = np.max(np.abs(spectrum))
        normalized = spectrum / scale if scale else spectrum
        options = dict(label=name, linestyle=style, marker=marker,
                       markevery=(5 * k, 28), markersize=4,
                       markerfacecolor="white", linewidth=1.1)
        axes[0].plot(omega, normalized.real, **options)
        axes[1].plot(omega, normalized.imag, **options)
    for ax, component in zip(axes, ["Real", "Imaginary"]):
        ax.set(xlabel=r"$\omega$ (inverse model-time units)",
               ylabel=f"{component} response / max |complex response|")
        ax.axhline(0, color="0.7", lw=0.6)
        ax.legend(fontsize=7)
    axes[0].set_title("Real part: all four models")
    axes[1].set_title("Imaginary part: all four models")
    print("Linear spectra:", len(spectra), "models; shape:", omega.shape)
    print("Curves are normalized separately; amplitudes are not compared.")
    plt.show()


if __name__ == "__main__":
    main()
