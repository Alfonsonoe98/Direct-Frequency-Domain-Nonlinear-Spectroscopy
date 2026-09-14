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
# 2. Laboratory time and pulse sequence
# ============================================================

dt = 0.05

times = np.arange(
    0.0,
    80.0 + dt,
    dt,
)

pulse_times = (
    5.0,
    15.0,
    35.0,
)

# Change this to any field-sign sector you want.
signature = (
    +1,
    -1,
    +1,
)

tau1, tau2, tau3 = pulse_times

t1_delay = (
    tau2
    - tau1
)

T_wait = (
    tau3
    - tau2
)


# ============================================================
# 3. Name the selected sector
# ============================================================

if signature == (+1, -1, +1):
    sector_name = "NR"
elif signature == (-1, +1, +1):
    sector_name = "R"
else:
    sector_name = str(signature)


# ============================================================
# 4. Continuous selected-pathway polarization
# ============================================================

P_pathway, components = (
    pulsus.impulsive_pathway_polarization(
        system=system,
        times=times,
        pulse_times=pulse_times,
        signature=signature,
        rho0=rho_ss,
        return_components=True,
    )
)

P1 = components["first_order"]
P2 = components["second_order"]
P3 = components["third_order"]


# ============================================================
# 5. Schematic Gaussian laser pulses
#
# IMPORTANT:
#
# The spectroscopy calculation above is IMPULSIVE.
# The true mathematical pulses are delta functions.
#
# These Gaussian fields are plotted only to visualize pulse
# positions, carrier oscillations, and an illustrative width.
# ============================================================

pulse_sigma = 1.0
omega_L = 5.00
pulse_amplitude = 1.0


def gaussian_pulse_field(
    times,
    center,
    sigma,
    omega_L,
    amplitude=1.0,
):
    delay = (
        times
        - center
    )

    envelope = (
        amplitude
        * np.exp(
            -0.5
            * (
                delay / sigma
            ) ** 2
        )
    )

    field = (
        envelope
        * np.cos(
            omega_L
            * delay
        )
    )

    return field, envelope


E1, envelope1 = gaussian_pulse_field(
    times,
    tau1,
    pulse_sigma,
    omega_L,
    pulse_amplitude,
)

E2, envelope2 = gaussian_pulse_field(
    times,
    tau2,
    pulse_sigma,
    omega_L,
    pulse_amplitude,
)

E3, envelope3 = gaussian_pulse_field(
    times,
    tau3,
    pulse_sigma,
    omega_L,
    pulse_amplitude,
)

E_total = (
    E1
    + E2
    + E3
)


# ============================================================
# 6. Diagnostics
# ============================================================

print()
print("IMPULSIVE PATHWAY POLARIZATION")
print()

print(
    "Pulse times:",
    pulse_times,
)

print(
    "Signature:",
    signature,
)

print(
    "Sector:",
    sector_name,
)

print(
    "t1 delay:",
    t1_delay,
)

print(
    "Waiting time T:",
    T_wait,
)

print()

print(
    "Max |P1|:",
    np.max(np.abs(P1)),
)

print(
    "Max |P2|:",
    np.max(np.abs(P2)),
)

print(
    "Max |P3|:",
    np.max(np.abs(P3)),
)

print()

print(
    "Illustrative Gaussian pulse sigma:",
    pulse_sigma,
)

print(
    "Carrier frequency omega_L:",
    omega_L,
)


# ============================================================
# 7. Normalize pathway using one common scale
# ============================================================

scale = np.max(
    np.abs(P_pathway)
)

if scale == 0.0:
    scale = 1.0

P_pathway_n = (
    P_pathway
    / scale
)

P1_n = (
    P1
    / scale
)

P2_n = (
    P2
    / scale
)

P3_n = (
    P3
    / scale
)


# ============================================================
# 8. Plot
# ============================================================

fig, axes = plt.subplots(
    3,
    1,
    figsize=(12, 10),
    sharex=True,
)


# ------------------------------------------------------------
# Panel 1: schematic laser pulses
# ------------------------------------------------------------

axes[0].plot(
    times,
    E_total,
    label=r"$E_1(t)+E_2(t)+E_3(t)$",
)

axes[0].plot(
    times,
    envelope1,
    linestyle="--",
)

axes[0].plot(
    times,
    envelope2,
    linestyle="--",
)

axes[0].plot(
    times,
    envelope3,
    linestyle="--",
)

axes[0].plot(
    times,
    -envelope1,
    linestyle="--",
)

axes[0].plot(
    times,
    -envelope2,
    linestyle="--",
)

axes[0].plot(
    times,
    -envelope3,
    linestyle="--",
)

axes[0].set_title(
    (
        "Schematic Gaussian laser pulses "
        "(impulsive calculation)"
    )
)

axes[0].set_ylabel(
    r"$E(t)$"
)

axes[0].legend()


# ------------------------------------------------------------
# Panel 2: continuous pathway polarization
# ------------------------------------------------------------

axes[1].plot(
    times,
    P_pathway_n.real,
    label=r"$\mathrm{Re}\,P(t)$",
)

axes[1].plot(
    times,
    P_pathway_n.imag,
    linestyle="--",
    label=r"$\mathrm{Im}\,P(t)$",
)

axes[1].axhline(
    0.0,
    linewidth=0.8,
)

axes[1].set_title(
    (
        "Continuous selected-pathway polarization: "
        f"{sector_name}"
    )
)

axes[1].set_ylabel(
    "Normalized polarization"
)

axes[1].legend()


# ------------------------------------------------------------
# Panel 3: active perturbative stage
# ------------------------------------------------------------

axes[2].plot(
    times,
    np.abs(P1_n),
    label=r"$|P^{(1)}|$",
)

axes[2].plot(
    times,
    np.abs(P2_n),
    label=r"$|P^{(2)}|$",
)

axes[2].plot(
    times,
    np.abs(P3_n),
    label=r"$|P^{(3)}|$",
)

axes[2].set_title(
    "Active pathway stage"
)

axes[2].set_xlabel(
    r"Laboratory time $t$"
)

axes[2].set_ylabel(
    "Magnitude"
)

axes[2].legend()


# ============================================================
# 9. Pulse markers
# ============================================================

for ax in axes:

    for number, pulse_time in enumerate(
        pulse_times,
        start=1,
    ):

        ax.axvline(
            pulse_time,
            linestyle=":",
        )

        ax.text(
            pulse_time,
            1.02,
            rf"$\tau_{number}$",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
        )


# ============================================================
# 10. Interval labels
# ============================================================

axes[1].text(
    0.5 * (tau1 + tau2),
    0.90,
    "1st-order state",
    transform=axes[1].get_xaxis_transform(),
    ha="center",
)

axes[1].text(
    0.5 * (tau2 + tau3),
    0.90,
    "2nd-order state",
    transform=axes[1].get_xaxis_transform(),
    ha="center",
)

axes[1].text(
    0.5 * (tau3 + times[-1]),
    0.90,
    "3rd-order state",
    transform=axes[1].get_xaxis_transform(),
    ha="center",
)


fig.suptitle(
    (
        "Impulsive pathway in laboratory time\n"
        rf"$t_1={t1_delay:g}$, "
        rf"$T={T_wait:g}$, "
        f"signature {signature}"
    ),
    fontsize=15,
)

fig.tight_layout()

plt.show()
