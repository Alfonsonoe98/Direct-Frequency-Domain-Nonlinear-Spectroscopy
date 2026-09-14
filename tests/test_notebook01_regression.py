import numpy as np

from scipy.integrate import solve_ivp
from numpy.linalg import norm

import pulsus

from dimer_model import build_dimer


def _vec(rho):
    return np.asarray(
        rho,
        dtype=complex,
    ).reshape(
        -1,
        order="F",
    )


def _unvec(rho_vec, d):
    return np.asarray(
        rho_vec,
        dtype=complex,
    ).reshape(
        (d, d),
        order="F",
    )


def _dissipator_rho(L, rho):
    LdagL = L.conj().T @ L

    return (
        L @ rho @ L.conj().T
        - 0.5
        * (
            LdagL @ rho
            + rho @ LdagL
        )
    )


def test_notebook01_impulsive_linear_td_matches_pulsus_fd():
    """
    Independent Notebook 01 regression.

    The time-domain calculation is performed directly in Hilbert
    space using the Lindblad master equation and solve_ivp.

    It does not use the PULSUS Liouvillian or time-domain response
    implementation.

    The resulting numerical Fourier transform is compared with the
    public PULSUS direct frequency-domain linear spectrum.
    """

    model = build_dimer()

    H = model["H_S"]
    mu = model["mu"]
    rho0 = model["rho0"]
    hbar = model["hbar"]
    d = model["d"]

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    # --------------------------------------------------------
    # Independent first-order initial condition
    #
    # rho^(1)(0) = (i / hbar) [mu, rho0]
    #
    # This is constructed directly in Hilbert space rather
    # than using the PULSUS interaction superoperator.
    # --------------------------------------------------------

    rho1_0 = (
        1j
        / hbar
        * (
            mu @ rho0
            - rho0 @ mu
        )
    )

    rho1_0_vec = _vec(
        rho1_0
    )

    # --------------------------------------------------------
    # Independent Hilbert-space Lindblad equation
    # --------------------------------------------------------

    def lindblad_rhs(t, rho_vec):
        rho = _unvec(
            rho_vec,
            d,
        )

        drho = (
            -1j
            / hbar
            * (
                H @ rho
                - rho @ H
            )
        )

        for L in collapse_ops:
            drho += _dissipator_rho(
                L,
                rho,
            )

        return _vec(
            drho
        )

    # Same high-accuracy transform window used in the
    # Notebook 01 validation.
    eta = 0.02

    dt = 0.02
    t_max = 250.0

    times = np.arange(
        0.0,
        t_max + dt,
        dt,
    )

    solution = solve_ivp(
        lindblad_rhs,
        t_span=(
            0.0,
            t_max,
        ),
        y0=rho1_0_vec,
        t_eval=times,
        method="DOP853",
        rtol=1.0e-10,
        atol=1.0e-12,
    )

    assert solution.success

    # --------------------------------------------------------
    # Independent time-domain polarization
    #
    # R^(1)(t) = Tr[mu rho^(1)(t)]
    # --------------------------------------------------------

    response_t = np.empty(
        len(times),
        dtype=complex,
    )

    for i, rho_vec in enumerate(
        solution.y.T
    ):
        rho = _unvec(
            rho_vec,
            d,
        )

        response_t[i] = np.trace(
            mu @ rho
        )

    # --------------------------------------------------------
    # Numerical causal Fourier transform
    # --------------------------------------------------------

    omega = np.linspace(
        0.6,
        1.4,
        81,
    )

    spectrum_td = np.empty(
        len(omega),
        dtype=complex,
    )

    for i, w in enumerate(
        omega
    ):
        integrand = (
            np.exp(
                (
                    1j * w
                    - eta
                )
                * times
            )
            * response_t
        )

        spectrum_td[i] = np.trapezoid(
            integrand,
            times,
        )

    # --------------------------------------------------------
    # PULSUS direct frequency-domain calculation
    # --------------------------------------------------------

    system = pulsus.SpectroscopySystem(
        H=H,
        dipole=mu,
        collapse_ops=collapse_ops,
        rho0=rho0,
        hbar=hbar,
    )

    spectrum_fd = (
        pulsus.impulsive_linear_spectrum(
            system=system,
            omega=omega,
            eta=eta,
            rho0=rho0,
        )
    )

    # --------------------------------------------------------
    # Notebook 01 found an error of approximately 1.2e-5.
    # Allow modest numerical variation between SciPy versions.
    # --------------------------------------------------------

    relative_error = (
        norm(
            spectrum_td
            - spectrum_fd
        )
        / norm(
            spectrum_fd
        )
    )

    assert relative_error < 5.0e-5
