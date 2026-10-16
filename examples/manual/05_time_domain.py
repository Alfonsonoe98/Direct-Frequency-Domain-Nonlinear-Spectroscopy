"""Impulsive RWA response functions and the selected pathway's laboratory trace."""

import matplotlib.pyplot as plt
import numpy as np
import pulsus

from common import build_dimer


def main():
    system, rho_ss = build_dimer()
    times = np.linspace(0.0, 80.0, 401)
    linear = pulsus.impulsive_linear_rwa_time_signal(
        system, times, sign=1, rho0=rho_ss
    )
    pulse_times = (5.0, 15.0, 35.0)
    pathway, components = pulsus.impulsive_pathway_polarization(
        system, times, pulse_times, pathway="NR", rho0=rho_ss,
        return_components=True,
    )
    np.testing.assert_allclose(pathway, sum(components.values()))
    assert np.all(pathway[times < pulse_times[0]] == 0)
    third_order = pulsus.impulsive_third_order_polarization(
        system, times, pulse_times, pathway="NR", rho0=rho_ss
    )
    after_third = times >= pulse_times[2]
    assert np.all(third_order[~after_third] == 0)
    np.testing.assert_allclose(pathway[after_third], third_order[after_third])
    assert np.isfinite(linear).all() and np.isfinite(pathway).all()
    delays = np.linspace(0, 30, 31)
    response = pulsus.impulsive_third_order_rwa_time_signal(
        system, t1=delays, t3=delays, T=20.0, pathway="NR", rho0=rho_ss
    )
    assert response.shape == (31, 31) and np.isfinite(response).all()
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.5), constrained_layout=True)
    axes[0].plot(times, linear.real, label="Real")
    axes[0].plot(times, linear.imag, label="Imaginary")
    axes[0].set(title="Impulsive linear RWA response", ylabel="Response")
    axes[0].legend()
    axes[1].plot(times, pathway.real, label="Real")
    axes[1].plot(times, pathway.imag, label="Imaginary")
    for t in pulse_times:
        axes[1].axvline(t, color="0.5", ls="--", lw=0.8)
    axes[1].set(title="Piecewise selected NR pathway (not total polarization)",
                xlabel="Laboratory time (model units)", ylabel="Pathway trace")
    axes[1].legend()
    print("Linear time trace shape:", linear.shape)
    print("Third-order time grid shape:", response.shape)
    print("Piecewise/third-order agreement after pulse 3: passed")
    print("Pathway component names:", sorted(components))
    plt.show()


if __name__ == "__main__":
    main()
