import numpy as np
import matplotlib.pyplot as plt

import pulsus
from dimer_model import build_dimer


# ============================================================
# 1. Frozen Notebook 04 reference
# ============================================================

reference = np.load(
    "tests/reference/notebook04_finite_pulse_reference.npz"
)

omega1_NR = reference["omega1_NR"]
omega1_R = reference["omega1_R"]
omega3 = reference["omega3"]

sigma = reference["sigma"].item()
T = reference["T"].item()
eta = reference["eta"].item()

omega_L1 = reference["omega_L1"].item()
omega_L2 = reference["omega_L2"].item()
omega_L3 = reference["omega_L3"].item()

E0 = reference["E0"].item()

phase1 = reference["phase1"].item()
phase2 = reference["phase2"].item()
phase3 = reference["phase3"].item()

ODE_NR = reference["NR_finite_ODE"]
ODE_R = reference["R_finite_ODE"]


# ============================================================
# 2. Build spectroscopy system
# ============================================================

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
    hbar=model["hbar"],
    dipole_plus=mu_plus,
    dipole_minus=mu_minus,
)

rho_ss = system.stationary_state()


# ============================================================
# 3. Gaussian pulses
# ============================================================

pulse1 = pulsus.GaussianPulse(
    omega_L=omega_L1,
    sigma=sigma,
    E0=E0,
    phase=phase1,
)

pulse2 = pulsus.GaussianPulse(
    omega_L=omega_L2,
    sigma=sigma,
    E0=E0,
    phase=phase2,
)

pulse3 = pulsus.GaussianPulse(
    omega_L=omega_L3,
    sigma=sigma,
    E0=E0,
    phase=phase3,
)


# ============================================================
# 4. Production PULSUS finite-pulse calculation
#
# F_j uses the physical L.
# G uses finite eta.
# ============================================================

FD_NR = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_NR,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="NR",
    rho0=rho_ss,
)

FD_R = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_R,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="R",
    rho0=rho_ss,
)


# ============================================================
# 5. Short-pulse RWA
# ============================================================

SHORT_NR = pulsus.short_pulse_rwa_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_NR,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="NR",
    rho0=rho_ss,
)

SHORT_R = pulsus.short_pulse_rwa_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_R,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="R",
    rho0=rho_ss,
)


# ============================================================
# 6. Impulsive RWA
# ============================================================

IMP_NR = pulsus.impulsive_rwa_spectrum(
    system=system,
    omega1=omega1_NR,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="NR",
    rho0=rho_ss,
)

IMP_R = pulsus.impulsive_rwa_spectrum(
    system=system,
    omega1=omega1_R,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="R",
    rho0=rho_ss,
)


# ============================================================
# 7. Complex spectral-shape discrepancy
# ============================================================

def shape_error(reference_spectrum, test_spectrum):
    ref = reference_spectrum.ravel()
    test = test_spectrum.ravel()

    alpha = (
        np.vdot(test, ref)
        / np.vdot(test, test)
    )

    residual = (
        ref
        - alpha * test
    )

    error = (
        np.linalg.norm(residual)
        / np.linalg.norm(ref)
    )

    return alpha, error


methods_NR = {
    "Finite-width PULSUS": FD_NR,
    "Short-pulse RWA": SHORT_NR,
    "Impulsive RWA": IMP_NR,
}

methods_R = {
    "Finite-width PULSUS": FD_R,
    "Short-pulse RWA": SHORT_R,
    "Impulsive RWA": IMP_R,
}


print()
print(
    f"sigma = {sigma}, T = {T}, eta = {eta}"
)

print()
print("NR errors relative to explicit ODE")

for name, spectrum in methods_NR.items():
    alpha, error = shape_error(
        ODE_NR,
        spectrum,
    )

    print(
        f"{name:22s}: "
        f"{100.0 * error:.6f} %"
        f"   alpha={alpha}"
    )


print()
print("R errors relative to explicit ODE")

for name, spectrum in methods_R.items():
    alpha, error = shape_error(
        ODE_R,
        spectrum,
    )

    print(
        f"{name:22s}: "
        f"{100.0 * error:.6f} %"
        f"   alpha={alpha}"
    )


# ============================================================
# 8. Exciton frequencies for guide lines
# ============================================================

exciton_energies = np.linalg.eigvalsh(
    model["H_S"]
)

# Remove ground and double-excitation energies.
# For this dimer, the two one-exciton energies are the middle
# two eigenvalues.
exciton_1 = exciton_energies[1]
exciton_2 = exciton_energies[2]

print()
print(
    "One-exciton frequencies:",
    exciton_1,
    exciton_2,
)


# ============================================================
# 9. Plot helper
# ============================================================

def normalized_imaginary(spectrum):
    data = np.imag(spectrum)

    scale = np.max(
        np.abs(data)
    )

    if scale > 0.0:
        data = data / scale

    return data


levels = np.linspace(
    -1.0,
    1.0,
    31,
)

NR_spectra = [
    ODE_NR,
    FD_NR,
    SHORT_NR,
    IMP_NR,
]

R_spectra = [
    ODE_R,
    FD_R,
    SHORT_R,
    IMP_R,
]

titles = [
    "Explicit finite-pulse ODE",
    "Finite-width PULSUS",
    "Short-pulse RWA",
    "Impulsive RWA",
]


# ============================================================
# 10. Figure
# ============================================================

fig, axes = plt.subplots(
    2,
    4,
    figsize=(14, 8),
    sharey="row",
)

last_contour = None

for col in range(4):

    # --------------------------------------------------------
    # NR row
    # --------------------------------------------------------

    data_NR = normalized_imaginary(
        NR_spectra[col]
    )

    last_contour = axes[0, col].contourf(
        omega1_NR,
        omega3,
        data_NR,
        levels=levels,
        cmap="RdBu_r",
    )

    axes[0, col].set_title(
        titles[col]
    )

    axes[0, col].set_xlabel(
        r"$\omega_1$"
    )

    axes[0, col].axvline(
        exciton_1,
        linestyle=":",
    )

    axes[0, col].axvline(
        exciton_2,
        linestyle=":",
    )

    axes[0, col].axhline(
        exciton_1,
        linestyle=":",
    )

    axes[0, col].axhline(
        exciton_2,
        linestyle=":",
    )

    # --------------------------------------------------------
    # R row
    # --------------------------------------------------------

    data_R = normalized_imaginary(
        R_spectra[col]
    )

    axes[1, col].contourf(
        omega1_R,
        omega3,
        data_R,
        levels=levels,
        cmap="RdBu_r",
    )

    axes[1, col].set_xlabel(
        r"$\nu_1$"
    )

    axes[1, col].axvline(
        -exciton_1,
        linestyle=":",
    )

    axes[1, col].axvline(
        -exciton_2,
        linestyle=":",
    )

    axes[1, col].axhline(
        exciton_1,
        linestyle=":",
    )

    axes[1, col].axhline(
        exciton_2,
        linestyle=":",
    )


axes[0, 0].set_ylabel(
    r"$\omega_3$"
)

axes[1, 0].set_ylabel(
    r"$\omega_3$"
)

fig.text(
    0.5,
    0.98,
    rf"NR comparison, $\sigma={sigma}$, $T={T}$",
    ha="center",
    va="top",
)

fig.text(
    0.5,
    0.49,
    rf"R comparison, $\sigma={sigma}$, $T={T}$",
    ha="center",
    va="top",
)

fig.colorbar(
    last_contour,
    ax=axes,
    label="Normalized imaginary response",
    shrink=0.75,
)

plt.show()
