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
# 4. Independent positive-detection sectors
# ============================================================

signatures = [
    (+1, +1, +1),
    (+1, -1, +1),   # NR
    (-1, +1, +1),   # R
    (-1, -1, +1),
]


results = {}


# ============================================================
# 5. Calculate full finite-pulse and short-pulse RWA
# ============================================================

for signature in signatures:

    s1, s2, s3 = signature

    omega1 = signed_axis(
        omega1_positive,
        s1,
    )

    full = pulsus.finite_pulse_spectrum(
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

    rwa = pulsus.short_pulse_rwa_spectrum(
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

    results[signature] = {
        "omega1": omega1,
        "full": full,
        "rwa": rwa,
    }


# ============================================================
# 6. Quantitative comparison
# ============================================================

global_full_max = max(
    np.max(np.abs(results[s]["full"]))
    for s in signatures
)


print()
print("FULL FINITE-PULSE vs SHORT-PULSE RWA")
print()


for signature in signatures:

    full = results[signature]["full"]
    rwa = results[signature]["rwa"]

    full_strength = (
        np.max(np.abs(full))
        / global_full_max
    )

    rwa_strength = (
        np.max(np.abs(rwa))
        / global_full_max
    )

    relative_error = (
        np.linalg.norm(full - rwa)
        / np.linalg.norm(full)
    )

    label = ""

    if signature == (+1, -1, +1):
        label = " NR"

    if signature == (-1, +1, +1):
        label = " R"

    print(
        f"{signature}{label}"
    )

    print(
        f"  Full strength: {full_strength:.6f}"
    )

    print(
        f"  RWA strength:  {rwa_strength:.6f}"
    )

    print(
        f"  Relative complex difference: "
        f"{100 * relative_error:.3f} %"
    )

    print()


# ============================================================
# 7. Common color scale
# ============================================================

global_imag_max = max(
    max(
        np.max(
            np.abs(
                np.imag(results[s]["full"])
            )
        ),
        np.max(
            np.abs(
                np.imag(results[s]["rwa"])
            )
        ),
    )
    for s in signatures
)

levels = np.linspace(
    -global_imag_max,
    global_imag_max,
    31,
)


# ============================================================
# 8. Plot
# ============================================================

fig, axes = plt.subplots(
    2,
    4,
    figsize=(16, 8),
)

contour = None


for column, signature in enumerate(signatures):

    omega1 = results[signature]["omega1"]
    full = results[signature]["full"]
    rwa = results[signature]["rwa"]

    label = str(signature)

    if signature == (+1, -1, +1):
        label += " [NR]"

    if signature == (-1, +1, +1):
        label += " [R]"

    # --------------------------------------------------------
    # Full finite-pulse
    # --------------------------------------------------------

    contour = axes[0, column].contourf(
        omega1,
        omega3,
        np.imag(full),
        levels=levels,
        cmap="RdBu_r",
        extend="both",
    )

    axes[0, column].set_title(
        label
    )

    axes[0, column].set_xlabel(
        r"$\omega_1$"
    )

    axes[0, column].set_ylabel(
        r"$\omega_3$"
    )

    # --------------------------------------------------------
    # Short-pulse RWA
    # --------------------------------------------------------

    axes[1, column].contourf(
        omega1,
        omega3,
        np.imag(rwa),
        levels=levels,
        cmap="RdBu_r",
        extend="both",
    )

    axes[1, column].set_xlabel(
        r"$\omega_1$"
    )

    axes[1, column].set_ylabel(
        r"$\omega_3$"
    )


axes[0, 0].text(
    -0.28,
    0.5,
    "Finite-pulse\nfull V",
    transform=axes[0, 0].transAxes,
    rotation=90,
    va="center",
    ha="center",
    fontsize=13,
)

axes[1, 0].text(
    -0.28,
    0.5,
    "Short-pulse\nRWA",
    transform=axes[1, 0].transAxes,
    rotation=90,
    va="center",
    ha="center",
    fontsize=13,
)


fig.colorbar(
    contour,
    ax=axes.ravel().tolist(),
    label="Imaginary response",
    shrink=0.85,
)

fig.suptitle(
    "Positive-detection third-order sectors",
    fontsize=15,
)

plt.show()
