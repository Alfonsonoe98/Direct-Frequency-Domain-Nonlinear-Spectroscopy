import numpy as np
import matplotlib.pyplot as plt

import pulsus
from dimer_model import build_dimer


# ============================================================
# 1. System
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
# 2. Pulses
# ============================================================

pulse1 = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
)

pulse2 = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
)

pulse3 = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
)

T = 20.0
eta = 0.02


# ============================================================
# 3. Frequency axes
# ============================================================

omega1_positive = np.linspace(
    0.7,
    1.5,
    81,
)

omega3 = np.linspace(
    0.7,
    1.5,
    81,
)


def signed_axis(axis_positive, sign):
    if sign == 1:
        return axis_positive.copy()

    return -axis_positive[::-1]


# ============================================================
# 4. Positive-detection sectors
# ============================================================

signatures = [
    (+1, +1, +1),
    (+1, -1, +1),   # NR
    (-1, +1, +1),   # R
    (-1, -1, +1),
]


# ============================================================
# 5. Helpers
# ============================================================

def finite_rwa_grid(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    signature,
    rho0,
):
    return pulsus.response_grid(
        response_function=pulsus.finite_pulse_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        T=T,
        eta=eta,
        signature=signature,
        rho0=rho0,
    )


def relative_error(reference, test):
    norm = np.linalg.norm(reference)

    if norm == 0.0:
        return 0.0

    return (
        np.linalg.norm(reference - test)
        / norm
    )


def shape_match(reference, test):
    """
    Find the optimal global complex factor alpha such that

        reference ~= alpha * test

    and return alpha and the residual shape error.
    """
    denominator = np.vdot(
        test,
        test,
    )

    if np.abs(denominator) == 0.0:
        return 0.0j, 1.0

    alpha = (
        np.vdot(test, reference)
        / denominator
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


# ============================================================
# 6. Calculate spectra
# ============================================================

results = {}


for signature in signatures:

    s1, s2, s3 = signature

    omega1 = signed_axis(
        omega1_positive,
        s1,
    )

    finite_full = pulsus.finite_pulse_spectrum(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        signature=signature,
        rho0=rho_ss,
    )

    finite_rwa = finite_rwa_grid(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        signature=signature,
        rho0=rho_ss,
    )

    short_rwa = pulsus.short_pulse_rwa_spectrum(
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        signature=signature,
        rho0=rho_ss,
    )

    impulsive_rwa = pulsus.impulsive_rwa_spectrum(
        system=system,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        signature=signature,
        rho0=rho_ss,
    )

    alpha_impulsive, impulsive_shape_error = shape_match(
        short_rwa,
        impulsive_rwa,
    )

    impulsive_scaled = (
        alpha_impulsive
        * impulsive_rwa
    )

    results[signature] = {
        "omega1": omega1,
        "finite_full": finite_full,
        "finite_rwa": finite_rwa,
        "short_rwa": short_rwa,
        "impulsive_rwa": impulsive_rwa,
        "impulsive_scaled": impulsive_scaled,
        "alpha_impulsive": alpha_impulsive,
        "impulsive_shape_error": impulsive_shape_error,
    }


# ============================================================
# 7. Quantitative comparison
# ============================================================

global_full_max = max(
    np.max(
        np.abs(
            results[s]["finite_full"]
        )
    )
    for s in signatures
)


print()
print("THIRD-ORDER APPROXIMATION LADDER")
print()


for signature in signatures:

    data = results[signature]

    finite_full = data["finite_full"]
    finite_rwa = data["finite_rwa"]
    short_rwa = data["short_rwa"]

    full_strength = (
        np.max(np.abs(finite_full))
        / global_full_max
    )

    finite_rwa_strength = (
        np.max(np.abs(finite_rwa))
        / global_full_max
    )

    short_strength = (
        np.max(np.abs(short_rwa))
        / global_full_max
    )

    rwa_error = relative_error(
        finite_full,
        finite_rwa,
    )

    short_error = relative_error(
        finite_rwa,
        short_rwa,
    )

    alpha_impulsive = data["alpha_impulsive"]
    impulsive_shape_error = data[
        "impulsive_shape_error"
    ]

    label = ""

    if signature == (+1, -1, +1):
        label = " NR"

    if signature == (-1, +1, +1):
        label = " R"

    print(
        f"{signature}{label}"
    )

    print(
        f"  Finite full V strength:   "
        f"{full_strength:.6f}"
    )

    print(
        f"  Finite RWA strength:      "
        f"{finite_rwa_strength:.6f}"
    )

    print(
        f"  Short-pulse RWA strength: "
        f"{short_strength:.6f}"
    )

    print(
        f"  Full -> finite RWA:       "
        f"{100 * rwa_error:.3f} %"
    )

    print(
        f"  Finite RWA -> short RWA:  "
        f"{100 * short_error:.3f} %"
    )

    print(
        f"  Short RWA -> impulsive "
        f"shape error: "
        f"{100 * impulsive_shape_error:.3f} %"
    )

    print(
        f"  Impulsive scale alpha:    "
        f"{alpha_impulsive}"
    )

    print()


# ============================================================
# 8. Plotting
#
# Rows 1-3 use the same absolute physical scale.
#
# Row 4 uses optimally scaled impulsive RWA because the
# impulsive molecular response does not contain the Gaussian
# field prefactors.
# ============================================================

plot_names = [
    "finite_full",
    "finite_rwa",
    "short_rwa",
    "impulsive_scaled",
]

row_titles = [
    "Finite-pulse full V",
    "Finite-pulse RWA",
    "Short-pulse RWA",
    "Impulsive RWA\n(shape-scaled)",
]


global_imag_max = max(
    np.max(
        np.abs(
            np.imag(
                results[s][name]
            )
        )
    )
    for s in signatures
    for name in plot_names
)

levels = np.linspace(
    -global_imag_max,
    global_imag_max,
    31,
)


fig, axes = plt.subplots(
    4,
    4,
    figsize=(16, 14),
)

contour = None


for column, signature in enumerate(signatures):

    omega1 = results[signature]["omega1"]

    label = str(signature)

    if signature == (+1, -1, +1):
        label += " [NR]"

    if signature == (-1, +1, +1):
        label += " [R]"

    axes[0, column].set_title(
        label
    )

    for row, name in enumerate(plot_names):

        spectrum = results[signature][name]

        contour = axes[row, column].contourf(
            omega1,
            omega3,
            np.imag(spectrum),
            levels=levels,
            cmap="RdBu_r",
            extend="both",
        )

        axes[row, column].set_xlabel(
            r"$\omega_1$"
        )

        axes[row, column].set_ylabel(
            r"$\omega_3$"
        )


for row, title in enumerate(row_titles):

    axes[row, 0].text(
        -0.34,
        0.5,
        title,
        transform=axes[row, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=12,
    )


fig.colorbar(
    contour,
    ax=axes.ravel().tolist(),
    label="Imaginary response",
    shrink=0.75,
)

fig.suptitle(
    "Third-order approximation ladder",
    fontsize=16,
)

plt.show()
