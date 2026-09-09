from pathlib import Path

import numpy as np

import pulsus
from dimer_model import build_dimer


REFERENCE_FILE = (
    Path(__file__).parent
    / "reference"
    / "notebook04_finite_pulse_reference.npz"
)


def finite_eta_response(
    L,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma,
    pathway,
    E0,
    phase1,
    phase2,
    phase3,
):
    """
    Reproduce the finite-eta pulse convention used only in
    Notebook 04's damped-ODE regression benchmark.
    """
    n = L.shape[0]
    I = np.eye(n, dtype=complex)

    # Notebook 04 uses L - eta I inside the pulse operators
    # for the consistent damped-transform comparison.
    L_pulse = L - eta * I

    s1, s2, s3 = pulsus.pathway_signs(pathway)

    F1 = pulsus.pulse_dressing_1(
        L_pulse,
        omega1=omega1,
        omega_L=omega_L1,
        sigma=sigma,
        E0=E0,
        phase=phase1,
        sign=s1,
    )

    F2 = pulsus.pulse_dressing_2(
        L_pulse,
        omega1=omega1,
        omega_L=omega_L2,
        sigma=sigma,
        E0=E0,
        phase=phase2,
        sign=s2,
    )

    F3 = pulsus.pulse_dressing_3(
        L_pulse,
        omega3=omega3,
        omega_L=omega_L3,
        sigma=sigma,
        E0=E0,
        phase=phase3,
        sign=s3,
    )

    return pulsus.third_order_response(
        L_super=L,
        V=V,
        observable=observable,
        rho0=rho0,
        omega1=omega1,
        omega3=omega3,
        T=T,
        eta=eta,
        F1=F1,
        F2=F2,
        F3=F3,
    )


def test_notebook04_finite_pulse_reference():
    reference = np.load(REFERENCE_FILE)

    model = build_dimer()

    L = model["L_super"]
    V = model["V"]
    mu = model["mu"]

    rho_ss = pulsus.stationary_state(
        L,
        model["d"],
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

    NR_reference = reference["NR_finite_FD"]
    R_reference = reference["R_finite_FD"]

    # Sample the 81 x 81 reference grid.
    # Full-grid comparison is unnecessary for every pytest run.
    indices = [0, 20, 40, 60, 80]

    for j in indices:
        for i in indices:

            NR = finite_eta_response(
                L=L,
                V=V,
                observable=mu,
                rho0=rho_ss,
                omega1=omega1_NR[i],
                omega3=omega3[j],
                T=T,
                eta=eta,
                omega_L1=omega_L1,
                omega_L2=omega_L2,
                omega_L3=omega_L3,
                sigma=sigma,
                pathway="NR",
                E0=E0,
                phase1=phase1,
                phase2=phase2,
                phase3=phase3,
            )

            R = finite_eta_response(
                L=L,
                V=V,
                observable=mu,
                rho0=rho_ss,
                omega1=omega1_R[i],
                omega3=omega3[j],
                T=T,
                eta=eta,
                omega_L1=omega_L1,
                omega_L2=omega_L2,
                omega_L3=omega_L3,
                sigma=sigma,
                pathway="R",
                E0=E0,
                phase1=phase1,
                phase2=phase2,
                phase3=phase3,
            )

            assert np.allclose(
                NR,
                NR_reference[j, i],
                rtol=1e-9,
                atol=1e-11,
            )

            assert np.allclose(
                R,
                R_reference[j, i],
                rtol=1e-9,
                atol=1e-11,
            )
