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

system = pulsus.SpectroscopySystem(
    H=model["H_S"],
    dipole=model["mu"],
    collapse_ops=collapse_ops,
    hbar=model["hbar"],
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
# 3. Positive-frequency base axes
# ============================================================

omega1_positive = np.linspace(
    0.7,
    1.5,
    81,
)

omega3_positive = np.linspace(
    0.7,
    1.5,
    81,
)


def signed_axis(axis_positive, sign):
    if sign == 1:
        return axis_positive.copy()

    return -axis_positive[::-1]


# ============================================================
# 4. All eight field-sign sectors
# ============================================================

signatures = [
    (+1, +1, +1),
    (+1, +1, -1),
    (+1, -1, +1),   # NR
    (+1, -1, -1),
    (-1, +1, +1),   # R
    (-1, +1, -1),
    (-1, -1, +1),
    (-1, -1, -1),
]

spectra = {}


for signature in signatures:

    s1, s2, s3 = signature

    omega1 = signed_axis(
        omega1_positive,
        s1,
    )

    omega3 = signed_axis(
        omega3_positive,
        s3,
    )

    spectrum = pulsus.finite_pulse_spectrum(
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

    spectra[signature] = {
        "omega1": omega1,
        "omega3": omega3,
        "spectrum": spectrum,
    }


# ============================================================
# 5. Relative strengths
# ============================================================

strengths = {
    signature: np.max(
        np.abs(
            spectra[signature]["spectrum"]
        )
    )
    for signature in signatures
}

global_strength = max(
    strengths.values()
)


print()
print("FINITE-PULSE FULL-V SIGNATURE STRENGTHS")
print()

for signature in signatures:

    relative = (
        strengths[signature]
        / global_strength
    )

    label = ""

    if signature == (+1, -1, +1):
        label = "  NR"

    if signature == (-1, +1, +1):
        label = "  R"

    print(
        f"{signature}: "
        f"{relative:.6f}"
        f"{label}"
    )


# ============================================================
# 6. Global plotting scale
#
# One common scale for all eight panels so that weak sectors
# actually look weak.
# ============================================================

global_imag_max = max(
    np.max(
        np.abs(
            np.imag(
                spectra[signature]["spectrum"]
            )
        )
    )
    for signature in signatures
)

levels = np.linspace(
    -global_imag_max,
    global_imag_max,
    31,
)


# ============================================================
# 7. Plot all eight sectors
# ============================================================

fig, axes = plt.subplots(
    2,
    4,
    figsize=(16, 8),
)

contour = None


for ax, signature in zip(
    axes.flat,
    signatures,
):

    omega1 = spectra[signature]["omega1"]
    omega3 = spectra[signature]["omega3"]
    spectrum = spectra[signature]["spectrum"]

    contour = ax.contourf(
        omega1,
        omega3,
        np.imag(spectrum),
        levels=levels,
        cmap="RdBu_r",
        extend="both",
    )

    relative = (
        strengths[signature]
        / global_strength
    )

    title = (
        f"{signature}\n"
        f"max |S| = {relative:.3f}"
    )

    if signature == (+1, -1, +1):
        title += "  [NR]"

    if signature == (-1, +1, +1):
        title += "  [R]"

    ax.set_title(
        title
    )

    ax.set_xlabel(
        r"$\omega_1$"
    )

    ax.set_ylabel(
        r"$\omega_3$"
    )


fig.colorbar(
    contour,
    ax=axes.ravel().tolist(),
    label="Imaginary response",
    shrink=0.85,
)

fig.suptitle(
    "Finite-pulse full-V third-order field-sign sectors",
    fontsize=15,
)

print()
print("SIGN-REVERSED PAIR RELATIONS")
print()

pairs = [
    ((+1, -1, +1), (-1, +1, -1)),
    ((-1, +1, +1), (+1, -1, -1)),
    ((+1, +1, +1), (-1, -1, -1)),
    ((+1, +1, -1), (-1, -1, +1)),
]

for signature, partner in pairs:

    S = spectra[signature]["spectrum"]
    S_partner = spectra[partner]["spectrum"]

    # Reverse both signed frequency axes so that
    # (-omega1, -omega3) aligns with (omega1, omega3).
    S_partner_mapped = S_partner[::-1, ::-1]

    error_plus = (
        np.linalg.norm(
            S_partner_mapped - np.conj(S)
        )
        / np.linalg.norm(S)
    )

    error_minus = (
        np.linalg.norm(
            S_partner_mapped + np.conj(S)
        )
        / np.linalg.norm(S)
    )

    print(
        signature,
        "<->",
        partner,
    )

    print(
        "  partner = +conj(S):",
        error_plus,
    )

    print(
        "  partner = -conj(S):",
        error_minus,
    )

plt.show()
