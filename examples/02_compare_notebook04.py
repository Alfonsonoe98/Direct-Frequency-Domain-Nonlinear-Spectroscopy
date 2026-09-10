import numpy as np

import pulsus
from dimer_model import build_dimer


# ============================================================
# 1. Load frozen Notebook 04 reference
# ============================================================

reference = np.load(
    "tests/reference/notebook04_finite_pulse_reference.npz"
)

omega1_NR = reference["omega1_NR"]
omega1_R = reference["omega1_R"]
omega3 = reference["omega3"]

sigma = reference["sigma"].item()
T = reference["T"].item()
eta = reference["eta"].item()

omega_L1 = reference["omega_L1"].item()
omega_L2 = reference["omega_L2"].item()
omega_L3 = reference["omega_L3"].item()

E0 = reference["E0"].item()

phase1 = reference["phase1"].item()
phase2 = reference["phase2"].item()
phase3 = reference["phase3"].item()


# ============================================================
# 2. Build the same dimer through PULSUS
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
# 3. Pulses
# ============================================================

pulse1 = pulsus.GaussianPulse(
    omega_L=omega_L1,
    sigma=sigma,
    E0=E0,
    phase=phase1,
)

pulse2 = pulsus.GaussianPulse(
    omega_L=omega_L2,
    sigma=sigma,
    E0=E0,
    phase=phase2,
)

pulse3 = pulsus.GaussianPulse(
    omega_L=omega_L3,
    sigma=sigma,
    E0=E0,
    phase=phase3,
)


# ============================================================
# 4. Production PULSUS finite-pulse spectra
#
# IMPORTANT:
# F_j uses L, with no finite-eta shift.
# Resolvents still use eta.
# ============================================================

S_NR = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_NR,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="NR",
    rho0=rho_ss,
)

S_R = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse1,
    pulse2=pulse2,
    pulse3=pulse3,
    omega1=omega1_R,
    omega3=omega3,
    T=T,
    eta=eta,
    pathway="R",
    rho0=rho_ss,
)


# ============================================================
# 5. Frozen explicit ODE references
# ============================================================

ODE_NR = reference["NR_finite_ODE"]
ODE_R = reference["R_finite_ODE"]


# ============================================================
# 6. Complex spectral-shape discrepancy
#
# Allow one optimal complex scalar alpha so that we compare
# spectral SHAPE rather than arbitrary global amplitude/phase.
# ============================================================

def shape_error(reference_spectrum, test_spectrum):
    ref = reference_spectrum.ravel()
    test = test_spectrum.ravel()

    alpha = (
        np.vdot(test, ref)
        / np.vdot(test, test)
    )

    residual = (
        ref
        - alpha * test
    )

    error = (
        np.linalg.norm(residual)
        / np.linalg.norm(ref)
    )

    return alpha, error


alpha_NR, error_NR = shape_error(
    ODE_NR,
    S_NR,
)

alpha_R, error_R = shape_error(
    ODE_R,
    S_R,
)


print("sigma =", sigma)
print("T     =", T)
print("eta   =", eta)

print()

print(
    "NR alpha:",
    alpha_NR,
)

print(
    "NR complex shape error:",
    100.0 * error_NR,
    "%",
)

print()

print(
    "R alpha:",
    alpha_R,
)

print(
    "R complex shape error:",
    100.0 * error_R,
    "%",
)
