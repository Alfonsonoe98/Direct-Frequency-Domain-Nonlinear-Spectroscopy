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
# 2. Time grids
#
# Use the same grid for t1 and t3.
# ============================================================

dt = 0.20
n_time = 2048

t1 = (
    np.arange(n_time)
    * dt
)

t3 = (
    np.arange(n_time)
    * dt
)

T = 20.0
eta = 0.02

signature = (
    +1,
    -1,
    +1,
)


# ============================================================
# 3. Time-domain NR response
#
# Shape:
#
#     (len(t3), len(t1))
#
# ============================================================

R_td = pulsus.impulsive_third_order_rwa_time_signal(
    system=system,
    t1=t1,
    t3=t3,
    T=T,
    signature=signature,
    rho0=rho_ss,
)


# ============================================================
# 4. Apply resolvent damping
#
# Direct PULSUS uses
#
# G_eta(omega)
#     = [(eta - i omega) I - L]^-1
#
# so the equivalent time-domain transform contains
#
# exp[-eta (t1 + t3)].
# ============================================================

damping = np.exp(
    -eta
    * (
        t3[:, np.newaxis]
        + t1[np.newaxis, :]
    )
)

R_damped = (
    damping
    * R_td
)


# ============================================================
# 5. Trapezoidal endpoint weights
# ============================================================

weighted = R_damped.copy()

weighted[0, :] *= 0.5
weighted[-1, :] *= 0.5

weighted[:, 0] *= 0.5
weighted[:, -1] *= 0.5


# ============================================================
# 6. Two-dimensional Fourier transform
#
# PULSUS convention:
#
# exp(+i omega1 t1)
# exp(+i omega3 t3)
#
# Therefore NumPy IFFT is the appropriate transform.
# ============================================================

S_fft_all = (
    dt
    * dt
    * n_time
    * n_time
    * np.fft.ifft2(
        weighted
    )
)


omega1_all = (
    2.0
    * np.pi
    * np.fft.fftfreq(
        n_time,
        d=dt,
    )
)

omega3_all = (
    2.0
    * np.pi
    * np.fft.fftfreq(
        n_time,
        d=dt,
    )
)


# ============================================================
# 7. Positive NR spectroscopy window
# ============================================================

mask1 = (
    (omega1_all >= 0.7)
    & (omega1_all <= 1.5)
)

mask3 = (
    (omega3_all >= 0.7)
    & (omega3_all <= 1.5)
)

omega1 = omega1_all[mask1]
omega3 = omega3_all[mask3]

S_fft = S_fft_all[
    np.ix_(
        mask3,
        mask1,
    )
]


# ============================================================
# 8. Direct frequency-domain NR spectrum
# ============================================================

S_fd = pulsus.impulsive_rwa_spectrum(
    system=system,
    omega1=omega1,
    omega3=omega3,
    T=T,
    eta=eta,
    signature=signature,
    rho0=rho_ss,
)


# ============================================================
# 9. Quantitative comparison
# ============================================================

relative_error = (
    np.linalg.norm(
        S_fd - S_fft
    )
    / np.linalg.norm(S_fd)
)


alpha = (
    np.vdot(
        S_fft,
        S_fd,
    )
    / np.vdot(
        S_fft,
        S_fft,
    )
)


shape_error = (
    np.linalg.norm(
        S_fd
        - alpha * S_fft
    )
    / np.linalg.norm(S_fd)
)


print()
print(
    "THIRD-ORDER NR IMPULSIVE RWA: "
    "TD+2D FFT vs DIRECT FD"
)
print()

print(
    "Time grid:",
    n_time,
    "x",
    n_time,
)

print(
    "dt:",
    dt,
)

print(
    "t_max:",
    t1[-1],
)

print(
    "T:",
    T,
)

print(
    "eta:",
    eta,
)

print(
    "FD grid:",
    len(omega3),
    "x",
    len(omega1),
)

print()

print(
    "Relative complex error:",
    100.0 * relative_error,
    "%",
)

print(
    "Complex shape error:",
    100.0 * shape_error,
    "%",
)

print(
    "Optimal alpha:",
    alpha,
)


# ============================================================
# 10. Plot
#
# Panel 1: time-domain response
# Panel 2: TD + FFT
# Panel 3: direct frequency domain
# ============================================================

fig, axes = plt.subplots(
    1,
    3,
    figsize=(15, 4.5),
)


# ------------------------------------------------------------
# Time domain
#
# Only show the early-time region so the structure is visible.
# ------------------------------------------------------------

time_plot_max = 40.0

time_mask1 = (
    t1 <= time_plot_max
)

time_mask3 = (
    t3 <= time_plot_max
)

R_plot = R_td[
    np.ix_(
        time_mask3,
        time_mask1,
    )
]


pulsus.plot_spectrum(
    t1[time_mask1],
    t3[time_mask3],
    R_plot,
    component="imag",
    normalize=True,
    ax=axes[0],
    colorbar=False,
    title="Time-domain NR response",
    xlabel=r"$t_1$",
    ylabel=r"$t_3$",
)


# ------------------------------------------------------------
# TD + 2D FFT
# ------------------------------------------------------------

pulsus.plot_spectrum(
    omega1,
    omega3,
    S_fft,
    component="imag",
    normalize=True,
    ax=axes[1],
    colorbar=False,
    title="TD + 2D FFT",
    xlabel=r"$\omega_1$",
    ylabel=r"$\omega_3$",
)


# ------------------------------------------------------------
# Direct frequency domain
# ------------------------------------------------------------

pulsus.plot_spectrum(
    omega1,
    omega3,
    S_fd,
    component="imag",
    normalize=True,
    ax=axes[2],
    colorbar=False,
    title="Direct FD",
    xlabel=r"$\omega_1$",
    ylabel=r"$\omega_3$",
)


fig.suptitle(
    "Impulsive NR RWA: time and frequency domains",
    fontsize=14,
)

fig.tight_layout()

plt.show()
