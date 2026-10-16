import matplotlib.pyplot as plt
import numpy as np
import pulsus

# Use consistent model units with hbar = 1.
# Build two coupled two-level sites.
identity = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)
a = np.kron(identity, lower)
b = np.kron(lower, identity)
H = (
    1.0 * (a.conj().T @ a)
    + 1.2 * (b.conj().T @ b)
    + 0.1 * (a.conj().T @ b + b.conj().T @ a)
)
dipole_plus = a.conj().T + b.conj().T
dipole_minus = dipole_plus.conj().T
dipole = dipole_plus + dipole_minus

gamma = 0.05
nbar = 0.05
collapse_ops = [
    np.sqrt(gamma * (nbar + 1)) * a,
    np.sqrt(gamma * nbar) * a.conj().T,
    np.sqrt(gamma * (nbar + 1)) * b,
    np.sqrt(gamma * nbar) * b.conj().T,
]
system = pulsus.SpectroscopySystem(
    H=H,
    dipole=dipole,
    collapse_ops=collapse_ops,
    hbar=1.0,
    dipole_plus=dipole_plus,
    dipole_minus=dipole_minus,
)
rho_ss = system.stationary_state()

# Define the pulse and frequency axes, then calculate the spectrum.
pulse = pulsus.GaussianPulse(omega_L=1.10, sigma=1.25, E0=1.0)
omega1 = np.linspace(0.7, 1.5, 31)
omega3 = np.linspace(0.7, 1.5, 31)
S = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    pathway="NR",
    rho0=rho_ss,
)

# Display the normalized imaginary component.
pulsus.plot_spectrum(
    omega1=omega1,
    omega3=omega3,
    spectrum=S,
    component="imag",
    normalize=True,
    title="Finite-pulse nonrephasing dimer spectrum",
    xlabel=r"$\omega_1$ (inverse time units)",
    ylabel=r"$\omega_3$ (inverse time units)",
)
plt.show()
