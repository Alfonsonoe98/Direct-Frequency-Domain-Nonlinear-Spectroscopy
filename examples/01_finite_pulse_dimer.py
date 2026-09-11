import numpy as np
import matplotlib.pyplot as plt

import pulsus


# ============================================================
# 1. Model parameters
# ============================================================

omega_a = 1.0
omega_b = 1.2
J = 0.1

gamma_a = 0.05
gamma_b = 0.05

nbar_a = 0.05
nbar_b = 0.05

mu_a = 1.0
mu_b = 1.0

hbar = 1.0


# ============================================================
# 2. Single-site operators
# ============================================================

I2 = np.eye(2, dtype=complex)

sigma_plus = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
    ],
    dtype=complex,
)

sigma_minus = sigma_plus.conj().T

# Tensor convention:
#
#     H_b tensor H_a
#
# giving basis
#
#     |gg>, |eg>, |ge>, |ee>
#

sigma_a_plus = np.kron(I2, sigma_plus)
sigma_a_minus = np.kron(I2, sigma_minus)

sigma_b_plus = np.kron(sigma_plus, I2)
sigma_b_minus = np.kron(sigma_minus, I2)


# ============================================================
# 3. Hamiltonian
# ============================================================

n_a = sigma_a_plus @ sigma_a_minus
n_b = sigma_b_plus @ sigma_b_minus

H = (
    omega_a * n_a
    + omega_b * n_b
    + J
    * (
        sigma_a_plus @ sigma_b_minus
        + sigma_b_plus @ sigma_a_minus
    )
)


# ============================================================
# 4. Dipole operator
# ============================================================

mu_plus = (
    mu_a * sigma_a_plus
    + mu_b * sigma_b_plus
)

mu_minus = (
    mu_a * sigma_a_minus
    + mu_b * sigma_b_minus
)

mu = mu_plus + mu_minus


# ============================================================
# 5. Lindblad collapse operators
# ============================================================

L_a_down = np.sqrt(
    gamma_a * (nbar_a + 1.0)
) * sigma_a_minus

L_a_up = np.sqrt(
    gamma_a * nbar_a
) * sigma_a_plus

L_b_down = np.sqrt(
    gamma_b * (nbar_b + 1.0)
) * sigma_b_minus

L_b_up = np.sqrt(
    gamma_b * nbar_b
) * sigma_b_plus

collapse_ops = [
    L_a_down,
    L_a_up,
    L_b_down,
    L_b_up,
]


# ============================================================
# 6. Build the spectroscopy system
# ============================================================

system = pulsus.SpectroscopySystem(
    H=H,
    dipole=mu,
    collapse_ops=collapse_ops,
    hbar=hbar,
    dipole_plus=mu_plus,
    dipole_minus=mu_minus,
)


# ============================================================
# 7. Stationary initial state
# ============================================================

rho_ss = system.stationary_state()

print(
    "Stationary-state trace:",
    np.trace(rho_ss),
)

print(
    "Stationary-state residual:",
    np.linalg.norm(
        system.L @ pulsus.vec(rho_ss)
    ),
)


# ============================================================
# 8. Gaussian laser pulses
# ============================================================

pulse = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
    E0=1.0,
    phase=0.0,
)


# ============================================================
# 9. Frequency grids
# ============================================================

omega1_NR = np.linspace(
    0.7,
    1.5,
    81,
)

# Native signed rephasing frequency axis.
#
# Keep it increasing from left to right:
#
#     -1.5 ... -0.7
#
omega1_R = -omega1_NR[::-1]

omega3 = np.linspace(
    0.7,
    1.5,
    81,
)

# ============================================================
# 10. Finite-pulse NR and R spectra
# ============================================================

S_NR = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1_NR,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    pathway="NR",
    rho0=rho_ss,
)

S_R = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1_R,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    pathway="R",
    rho0=rho_ss,
)

print(
    "NR spectrum shape:",
    S_NR.shape,
)

print(
    "R spectrum shape:",
    S_R.shape,
)

print(
    "Maximum |S_NR|:",
    np.max(np.abs(S_NR)),
)

print(
    "Maximum |S_R|:",
    np.max(np.abs(S_R)),
)

# ============================================================
# 11. Check selected points against the low-level API
# ============================================================

check_indices = [
    (0, 0),
    (20, 20),
    (40, 40),
    (60, 60),
    (80, 80),
    (20, 60),
    (60, 20),
]

max_difference_NR = 0.0
max_difference_R = 0.0

for j, i in check_indices:

    direct_NR = pulsus.finite_pulse_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=rho_ss,
        omega1=omega1_NR[i],
        omega3=omega3[j],
        T=20.0,
        eta=0.02,
        omega_L1=pulse.omega_L,
        omega_L2=pulse.omega_L,
        omega_L3=pulse.omega_L,
        sigma1=pulse.sigma,
        sigma2=pulse.sigma,
        sigma3=pulse.sigma,
        pathway="NR",
        E01=pulse.E0,
        E02=pulse.E0,
        E03=pulse.E0,
        phase1=pulse.phase,
        phase2=pulse.phase,
        phase3=pulse.phase,
    )

    direct_R = pulsus.finite_pulse_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=rho_ss,
        omega1=omega1_R[i],
        omega3=omega3[j],
        T=20.0,
        eta=0.02,
        omega_L1=pulse.omega_L,
        omega_L2=pulse.omega_L,
        omega_L3=pulse.omega_L,
        sigma1=pulse.sigma,
        sigma2=pulse.sigma,
        sigma3=pulse.sigma,
        pathway="R",
        E01=pulse.E0,
        E02=pulse.E0,
        E03=pulse.E0,
        phase1=pulse.phase,
        phase2=pulse.phase,
        phase3=pulse.phase,
    )

    max_difference_NR = max(
        max_difference_NR,
        abs(S_NR[j, i] - direct_NR),
    )

    max_difference_R = max(
        max_difference_R,
        abs(S_R[j, i] - direct_R),
    )

print(
    "Maximum NR public/low-level difference:",
    max_difference_NR,
)

print(
    "Maximum R public/low-level difference:",
    max_difference_R,
)

# ============================================================
# 13. Plot NR spectrum
# ============================================================

fig_NR, ax_NR, contour_NR = pulsus.plot_spectrum(
    omega1=omega1_NR,
    omega3=omega3,
    spectrum=S_NR,
    component="imag",
    normalize=True,
    title="Finite-pulse NR spectrum",
)

fig_NR.tight_layout()


# ============================================================
# 14. Plot R spectrum
# ============================================================

fig_R, ax_R, contour_R = pulsus.plot_spectrum(
    omega1=omega1_R,
    omega3=omega3,
    spectrum=S_R,
    component="imag",
    normalize=True,
    title="Finite-pulse R spectrum",
)

fig_R.tight_layout()


plt.show()
