import numpy as np
import matplotlib.pyplot as plt

import pulsus
from dimer_model import build_dimer


# ============================================================
# 1. Coupled-dimer system
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
# 2. Pulse and frequency grid
# ============================================================

pulse = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
    E0=1.0,
    phase=0.0,
)

omega = np.linspace(
    0.6,
    1.4,
    801,
)

eta = 0.02


# ============================================================
# 3. Impulsive full-V and RWA responses
# ============================================================

M_full = pulsus.impulsive_linear_spectrum(
    system=system,
    omega=omega,
    eta=eta,
    rho0=rho_ss,
)

M_rwa = pulsus.impulsive_linear_rwa_spectrum(
    system=system,
    omega=omega,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)


# ============================================================
# 4. Gaussian-driven full-V and RWA responses
# ============================================================

P_full = pulsus.finite_pulse_linear_spectrum(
    system=system,
    pulse=pulse,
    omega=omega,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)

P_rwa = pulsus.short_pulse_linear_rwa_spectrum(
    system=system,
    pulse=pulse,
    omega=omega,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)


# ============================================================
# 5. Quantitative errors
# ============================================================

def relative_error(reference, test):
    return (
        np.linalg.norm(reference - test)
        / np.linalg.norm(reference)
    )


def shape_error(reference, test):
    alpha = (
        np.vdot(test, reference)
        / np.vdot(test, test)
    )

    residual = (
        reference
        - alpha * test
    )

    error = (
        np.linalg.norm(residual)
        / np.linalg.norm(reference)
    )

    return alpha, error


alpha_M, error_M_shape = shape_error(
    M_full,
    M_rwa,
)

alpha_P, error_P_shape = shape_error(
    P_full,
    P_rwa,
)

error_M = relative_error(
    M_full,
    M_rwa,
)

error_P = relative_error(
    P_full,
    P_rwa,
)


print()
print("POSITIVE-FREQUENCY LINEAR RESPONSE")
print()

print("Impulsive full V vs RWA")
print(
    "  Relative complex error:",
    100.0 * error_M,
    "%",
)
print(
    "  Complex shape error:",
    100.0 * error_M_shape,
    "%",
)
print(
    "  Optimal alpha:",
    alpha_M,
)

print()

print("Gaussian full V vs RWA")
print(
    "  Relative complex error:",
    100.0 * error_P,
    "%",
)
print(
    "  Complex shape error:",
    100.0 * error_P_shape,
    "%",
)
print(
    "  Optimal alpha:",
    alpha_P,
)


# ============================================================
# 6. Exciton frequencies
# ============================================================

energies = np.linalg.eigvalsh(
    model["H_S"]
)

omega_minus = energies[1]
omega_plus = energies[2]


# ============================================================
# 7. Common normalization
#
# Full-V and RWA curves in each row use the same scale.
# ============================================================

scale_M = np.max(
    np.abs(M_full)
)

if scale_M == 0.0:
    scale_M = 1.0

scale_P = np.max(
    np.abs(P_full)
)

if scale_P == 0.0:
    scale_P = 1.0


M_full_n = M_full / scale_M
M_rwa_n = M_rwa / scale_M

P_full_n = P_full / scale_P
P_rwa_n = P_rwa / scale_P


delta_M = (
    np.abs(M_full - M_rwa)
    / scale_M
)

delta_P = (
    np.abs(P_full - P_rwa)
    / scale_P
)


# ============================================================
# 8. Figure
# ============================================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(15, 8),
    sharex=True,
)


# ------------------------------------------------------------
# Row 1: impulsive
# ------------------------------------------------------------

axes[0, 0].plot(
    omega,
    M_full_n.imag,
    label="Full V",
)

axes[0, 0].plot(
    omega,
    M_rwa_n.imag,
    linestyle="--",
    label="RWA",
)

axes[0, 0].set_title(
    "Impulsive: imaginary part"
)

axes[0, 0].set_ylabel(
    r"$\mathrm{Im}\,M^{(1)} / \max |M_{\rm full}^{(1)}|$"
)

axes[0, 0].legend()


axes[0, 1].plot(
    omega,
    M_full_n.real,
    label="Full V",
)

axes[0, 1].plot(
    omega,
    M_rwa_n.real,
    linestyle="--",
    label="RWA",
)

axes[0, 1].set_title(
    "Impulsive: real part"
)

axes[0, 1].set_ylabel(
    r"$\mathrm{Re}\,M^{(1)} / \max |M_{\rm full}^{(1)}|$"
)

axes[0, 1].legend()


axes[0, 2].plot(
    omega,
    delta_M,
)

axes[0, 2].set_title(
    "Impulsive: full V - RWA"
)

axes[0, 2].set_ylabel(
    r"$|M_{\rm full}^{(1)}-M_{\rm RWA}^{(1)}|"
    r"/\max|M_{\rm full}^{(1)}|$"
)


# ------------------------------------------------------------
# Row 2: Gaussian-driven
# ------------------------------------------------------------

axes[1, 0].plot(
    omega,
    P_full_n.imag,
    label="Full V",
)

axes[1, 0].plot(
    omega,
    P_rwa_n.imag,
    linestyle="--",
    label="RWA",
)

axes[1, 0].set_title(
    "Gaussian-driven: imaginary part"
)

axes[1, 0].set_ylabel(
    r"$\mathrm{Im}\,P^{(1)} / \max |P_{\rm full}^{(1)}|$"
)

axes[1, 0].legend()


axes[1, 1].plot(
    omega,
    P_full_n.real,
    label="Full V",
)

axes[1, 1].plot(
    omega,
    P_rwa_n.real,
    linestyle="--",
    label="RWA",
)

axes[1, 1].set_title(
    "Gaussian-driven: real part"
)

axes[1, 1].set_ylabel(
    r"$\mathrm{Re}\,P^{(1)} / \max |P_{\rm full}^{(1)}|$"
)

axes[1, 1].legend()


axes[1, 2].plot(
    omega,
    delta_P,
)

axes[1, 2].set_title(
    "Gaussian-driven: full V - RWA"
)

axes[1, 2].set_ylabel(
    r"$|P_{\rm full}^{(1)}-P_{\rm RWA}^{(1)}|"
    r"/\max|P_{\rm full}^{(1)}|$"
)


# ============================================================
# 9. Exciton guide lines and labels
# ============================================================

for ax in axes.flat:

    ax.axvline(
        omega_minus,
        linestyle=":",
    )

    ax.axvline(
        omega_plus,
        linestyle=":",
    )

    ax.set_xlabel(
        r"$\omega$"
    )


fig.suptitle(
    "Linear response: full interaction vs RWA",
    fontsize=15,
)

fig.tight_layout()

plt.show()
