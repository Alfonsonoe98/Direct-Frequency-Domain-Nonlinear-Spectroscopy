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
# 2. Time grid
#
# The direct frequency-domain resolvent uses
#
# G_eta(omega) = [(eta - i omega) I - L]^{-1}
#
# corresponding to
#
# integral_0^inf dt exp(-eta t) exp(+i omega t) R(t).
#
# NumPy's IFFT uses the required positive exponential.
# ============================================================

dt = 0.10
n_time = 4096

times = (
    np.arange(n_time)
    * dt
)

eta = 0.02


# ============================================================
# 3. Impulsive RWA response in time
# ============================================================

R_t = pulsus.impulsive_linear_rwa_time_signal(
    system=system,
    times=times,
    sign=+1,
    rho0=rho_ss,
)


# ============================================================
# 4. Apply the same eta damping used by the resolvent
# ============================================================

R_damped = (
    np.exp(-eta * times)
    * R_t
)


# ============================================================
# 5. FFT
#
# We use an IFFT because PULSUS uses the Fourier convention
#
#     exp(+i omega t).
#
# Zero padding gives a smoother frequency grid.
#
# Endpoint weights implement the trapezoidal rule.
# ============================================================

weighted = R_damped.copy()

weighted[0] *= 0.5
weighted[-1] *= 0.5

n_fft = 32768

M_fft_all = (
    dt
    * n_fft
    * np.fft.ifft(
        weighted,
        n=n_fft,
    )
)

omega_all = (
    2.0
    * np.pi
    * np.fft.fftfreq(
        n_fft,
        d=dt,
    )
)


# ============================================================
# 6. Keep the positive-frequency spectroscopy window
# ============================================================

mask = (
    (omega_all >= 0.6)
    & (omega_all <= 1.4)
)

omega = omega_all[mask]
M_fft = M_fft_all[mask]


# ============================================================
# 7. Direct PULSUS frequency-domain calculation
# ============================================================

M_fd = pulsus.impulsive_linear_rwa_spectrum(
    system=system,
    omega=omega,
    eta=eta,
    sign=+1,
    rho0=rho_ss,
)


# ============================================================
# 8. Quantitative comparison
# ============================================================

relative_error = (
    np.linalg.norm(
        M_fd - M_fft
    )
    / np.linalg.norm(M_fd)
)

alpha = (
    np.vdot(
        M_fft,
        M_fd,
    )
    / np.vdot(
        M_fft,
        M_fft,
    )
)

shape_error = (
    np.linalg.norm(
        M_fd
        - alpha * M_fft
    )
    / np.linalg.norm(M_fd)
)


print()
print("LINEAR IMPULSIVE RWA: TD+FFT vs DIRECT FD")
print()

print(
    "Number of time points:",
    n_time,
)

print(
    "dt:",
    dt,
)

print(
    "t_max:",
    times[-1],
)

print(
    "eta:",
    eta,
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
# 9. Plot
# ============================================================

scale = np.max(
    np.abs(M_fd)
)

M_fd_n = (
    M_fd
    / scale
)

M_fft_n = (
    M_fft
    / scale
)


fig, axes = plt.subplots(
    1,
    2,
    figsize=(11, 4.5),
)


# ------------------------------------------------------------
# Imaginary part
# ------------------------------------------------------------

axes[0].plot(
    omega,
    M_fd_n.imag,
    label="Direct FD",
)

axes[0].plot(
    omega,
    M_fft_n.imag,
    linestyle="--",
    label="TD + FFT",
)

axes[0].set_title(
    "Imaginary response"
)

axes[0].set_xlabel(
    r"$\omega$"
)

axes[0].set_ylabel(
    "Normalized response"
)

axes[0].legend()


# ------------------------------------------------------------
# Real part
# ------------------------------------------------------------

axes[1].plot(
    omega,
    M_fd_n.real,
    label="Direct FD",
)

axes[1].plot(
    omega,
    M_fft_n.real,
    linestyle="--",
    label="TD + FFT",
)

axes[1].set_title(
    "Real response"
)

axes[1].set_xlabel(
    r"$\omega$"
)

axes[1].set_ylabel(
    "Normalized response"
)

axes[1].legend()


fig.suptitle(
    "Impulsive linear RWA: time domain vs direct frequency domain",
    fontsize=14,
)

fig.tight_layout()

plt.show()
