import numpy as np
import matplotlib.pyplot as plt

import pulsus
from dimer_model import build_dimer


# ============================================================
# 1. Build the same coupled-dimer spectroscopy system
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
# 2. Gaussian pulse
# ============================================================

pulse = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
    E0=1.0,
    phase=0.0,
)


# ============================================================
# 3. Signed frequency axes
# ============================================================

omega_positive = np.linspace(
    0.6,
    1.4,
    401,
)

omega_negative = np.linspace(
    -1.4,
    -0.6,
    401,
)

eta = 0.02


# ============================================================
# 4. Positive-frequency sector
# ============================================================

# Impulsive, full molecular interaction
M_full_positive = pulsus.impulsive_linear_spectrum(
    system=system,
    omega=omega_positive,
    eta=eta,
    rho0=rho_ss,
)

# Impulsive RWA, positive-frequency molecular interaction
M_rwa_positive = pulsus.impulsive_linear_rwa_spectrum(
    system=system,
    omega=omega_positive,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)

# Gaussian field component with full molecular interaction
P_full_positive = pulsus.finite_pulse_linear_spectrum(
    system=system,
    pulse=pulse,
    omega=omega_positive,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)

# Gaussian field component with molecular RWA
P_rwa_positive = pulsus.short_pulse_linear_rwa_spectrum(
    system=system,
    pulse=pulse,
    omega=omega_positive,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)


# ============================================================
# 5. Negative-frequency sector
# ============================================================

M_full_negative = pulsus.impulsive_linear_spectrum(
    system=system,
    omega=omega_negative,
    eta=eta,
    rho0=rho_ss,
)

M_rwa_negative = pulsus.impulsive_linear_rwa_spectrum(
    system=system,
    omega=omega_negative,
    eta=eta,
    sign=-1,
    rho0=rho_ss,
)

P_full_negative = pulsus.finite_pulse_linear_spectrum(
    system=system,
    pulse=pulse,
    omega=omega_negative,
    eta=eta,
    sign=-1,
    rho0=rho_ss,
)

P_rwa_negative = pulsus.short_pulse_linear_rwa_spectrum(
    system=system,
    pulse=pulse,
    omega=omega_negative,
    eta=eta,
    sign=-1,
    rho0=rho_ss,
)


# ============================================================
# 6. Quantitative comparison
# ============================================================

def relative_error(reference, test):
    return (
        np.linalg.norm(reference - test)
        / np.linalg.norm(reference)
    )


def shape_error(reference, test):
    """
    Complex spectral-shape error after allowing one optimal
    global complex scale factor.
    """
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


comparisons = [
    (
        "Positive impulsive",
        M_full_positive,
        M_rwa_positive,
    ),
    (
        "Positive Gaussian",
        P_full_positive,
        P_rwa_positive,
    ),
    (
        "Negative impulsive",
        M_full_negative,
        M_rwa_negative,
    ),
    (
        "Negative Gaussian",
        P_full_negative,
        P_rwa_negative,
    ),
]


print()
print("FULL V vs RWA LINEAR RESPONSE")
print()

for name, full, rwa in comparisons:

    rel = relative_error(
        full,
        rwa,
    )

    alpha, shape = shape_error(
        full,
        rwa,
    )

    print(name)

    print(
        "  Relative complex error:",
        100.0 * rel,
        "%",
    )

    print(
        "  Complex shape error:",
        100.0 * shape,
        "%",
    )

    print(
        "  Optimal alpha:",
        alpha,
    )

    print()


# ============================================================
# 7. Exciton frequencies
# ============================================================

energies = np.linalg.eigvalsh(
    model["H_S"]
)

omega_minus = energies[1]
omega_plus = energies[2]

print(
    "One-exciton frequencies:",
    omega_minus,
    omega_plus,
)


# ============================================================
# 8. Plotting helper
#
# Both curves in each panel are normalized by the SAME
# full-V scale. RWA amplitude differences therefore remain
# visible.
# ============================================================

def normalized_imag_pair(
    reference,
    test,
):
    scale = np.max(
        np.abs(
            np.imag(reference)
        )
    )

    if scale == 0.0:
        scale = 1.0

    return (
        np.imag(reference) / scale,
        np.imag(test) / scale,
    )


Mi_full_pos, Mi_rwa_pos = normalized_imag_pair(
    M_full_positive,
    M_rwa_positive,
)

Pi_full_pos, Pi_rwa_pos = normalized_imag_pair(
    P_full_positive,
    P_rwa_positive,
)

Mi_full_neg, Mi_rwa_neg = normalized_imag_pair(
    M_full_negative,
    M_rwa_negative,
)

Pi_full_neg, Pi_rwa_neg = normalized_imag_pair(
    P_full_negative,
    P_rwa_negative,
)


# ============================================================
# 9. Comparison figure
# ============================================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(11, 8),
)

# ------------------------------------------------------------
# Positive impulsive
# ------------------------------------------------------------

axes[0, 0].plot(
    omega_positive,
    Mi_full_pos,
    label="Full V",
)

axes[0, 0].plot(
    omega_positive,
    Mi_rwa_pos,
    linestyle="--",
    label="RWA",
)

axes[0, 0].set_title(
    "Impulsive linear response, +ω sector"
)

axes[0, 0].set_xlabel(
    r"$\omega$"
)

axes[0, 0].set_ylabel(
    "Normalized imaginary response"
)

axes[0, 0].legend()


# ------------------------------------------------------------
# Positive Gaussian
# ------------------------------------------------------------

axes[0, 1].plot(
    omega_positive,
    Pi_full_pos,
    label="Full V + Gaussian field",
)

axes[0, 1].plot(
    omega_positive,
    Pi_rwa_pos,
    linestyle="--",
    label="RWA + Gaussian field",
)

axes[0, 1].set_title(
    "Gaussian-driven linear response, +ω sector"
)

axes[0, 1].set_xlabel(
    r"$\omega$"
)

axes[0, 1].set_ylabel(
    "Normalized imaginary response"
)

axes[0, 1].legend()


# ------------------------------------------------------------
# Negative impulsive
# ------------------------------------------------------------

axes[1, 0].plot(
    omega_negative,
    Mi_full_neg,
    label="Full V",
)

axes[1, 0].plot(
    omega_negative,
    Mi_rwa_neg,
    linestyle="--",
    label="RWA",
)

axes[1, 0].set_title(
    "Impulsive linear response, -ω sector"
)

axes[1, 0].set_xlabel(
    r"$\omega$"
)

axes[1, 0].set_ylabel(
    "Normalized imaginary response"
)

axes[1, 0].legend()


# ------------------------------------------------------------
# Negative Gaussian
# ------------------------------------------------------------

axes[1, 1].plot(
    omega_negative,
    Pi_full_neg,
    label="Full V + Gaussian field",
)

axes[1, 1].plot(
    omega_negative,
    Pi_rwa_neg,
    linestyle="--",
    label="RWA + Gaussian field",
)

axes[1, 1].set_title(
    "Gaussian-driven linear response, -ω sector"
)

axes[1, 1].set_xlabel(
    r"$\omega$"
)

axes[1, 1].set_ylabel(
    "Normalized imaginary response"
)

axes[1, 1].legend()


# ============================================================
# 10. Exciton guide lines
# ============================================================

for ax in axes[0, :]:
    ax.axvline(
        omega_minus,
        linestyle=":",
    )

    ax.axvline(
        omega_plus,
        linestyle=":",
    )


for ax in axes[1, :]:
    ax.axvline(
        -omega_plus,
        linestyle=":",
    )

    ax.axvline(
        -omega_minus,
        linestyle=":",
    )


fig.tight_layout()

plt.show()
