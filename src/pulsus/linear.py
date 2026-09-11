import numpy as np

from .liouville import vec
from .pulses import (
    gaussian_spectrum,
    pulse_dressing_1,
)
from .resolvents import resolvent_action
from .response import observable_bra


def _validate_sign(sign):
    if sign not in (-1, 1):
        raise ValueError(
            "sign must be +1 or -1"
        )

    return sign


def _rwa_interaction(
    V_plus,
    V_minus,
    sign,
):
    sign = _validate_sign(sign)

    if sign == 1:
        return V_plus

    return V_minus


# ============================================================
# Low-level full-interaction linear response
# ============================================================

def impulsive_linear_response(
    L_super,
    V,
    observable,
    rho0,
    omega,
    eta,
):
    """
    Impulsive first-order molecular response using the full
    light-matter interaction.

    M^(1)(omega)
        = <<observable| G_eta(omega) V |rho0>>.
    """
    rho_vec = vec(rho0)

    state = np.asarray(
        V,
        dtype=complex,
    ) @ rho_vec

    state = resolvent_action(
        L_super=L_super,
        omega=omega,
        rhs=state,
        eta=eta,
    )

    return (
        observable_bra(observable)
        @ state
    )


def finite_pulse_linear_response(
    L_super,
    V,
    observable,
    rho0,
    omega,
    eta,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    First-order response to one Gaussian field-frequency
    component while retaining the full molecular interaction V.

    P^(1,s)(omega)
        = <<observable|
          G_eta(omega)
          V
          F_1^(s)(omega)
          |rho0>>.

    This is NOT a molecular RWA calculation.  The sign selects
    the laser field component, while V remains complete.
    """
    sign = _validate_sign(sign)

    F1 = pulse_dressing_1(
        L_super=L_super,
        omega1=omega,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )

    state = F1 @ vec(rho0)
    state = np.asarray(
        V,
        dtype=complex,
    ) @ state

    state = resolvent_action(
        L_super=L_super,
        omega=omega,
        rhs=state,
        eta=eta,
    )

    return (
        observable_bra(observable)
        @ state
    )


# ============================================================
# Low-level RWA linear response
# ============================================================

def impulsive_linear_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega,
    eta,
    sign=1,
):
    """
    Impulsive first-order RWA response for a selected
    field-frequency sign.

    sign = +1 -> V_plus
    sign = -1 -> V_minus
    """
    V_s = _rwa_interaction(
        V_plus,
        V_minus,
        sign,
    )

    state = V_s @ vec(rho0)

    state = resolvent_action(
        L_super=L_super,
        omega=omega,
        rhs=state,
        eta=eta,
    )

    return (
        observable_bra(observable)
        @ state
    )


def short_pulse_linear_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega,
    eta,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    Short-pulse RWA first-order polarization.

    P_short^(1,s)(omega)
        = E^(s)(omega)
          M_RWA^(1,s)(omega).
    """
    sign = _validate_sign(sign)

    molecular = impulsive_linear_rwa_response(
        L_super=L_super,
        V_plus=V_plus,
        V_minus=V_minus,
        observable=observable,
        rho0=rho0,
        omega=omega,
        eta=eta,
        sign=sign,
    )

    field = gaussian_spectrum(
        omega=omega,
        omega_L=omega_L,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )

    return field * molecular


# ============================================================
# High-level signal API
# ============================================================

def impulsive_linear_signal(
    system,
    omega,
    eta,
    rho0=None,
):
    """
    High-level impulsive full-interaction linear response.
    """
    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    return impulsive_linear_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=rho0,
        omega=omega,
        eta=eta,
    )


def finite_pulse_linear_signal(
    system,
    pulse,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    High-level finite-Gaussian full-interaction response.

    The field sign is selected explicitly, but the complete
    molecular interaction V is retained.
    """
    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    return finite_pulse_linear_response(
        L_super=system.L,
        V=system.V,
        observable=system.dipole,
        rho0=rho0,
        omega=omega,
        eta=eta,
        omega_L=pulse.omega_L,
        sigma=pulse.sigma,
        E0=pulse.E0,
        phase=pulse.phase,
        sign=sign,
    )


def impulsive_linear_rwa_signal(
    system,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    High-level impulsive RWA linear response.
    """
    if (
        system.V_plus is None
        or system.V_minus is None
    ):
        raise ValueError(
            "RWA calculations require dipole_plus and "
            "dipole_minus in SpectroscopySystem"
        )

    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    return impulsive_linear_rwa_response(
        L_super=system.L,
        V_plus=system.V_plus,
        V_minus=system.V_minus,
        observable=system.dipole,
        rho0=rho0,
        omega=omega,
        eta=eta,
        sign=sign,
    )


def short_pulse_linear_rwa_signal(
    system,
    pulse,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    High-level short-pulse RWA linear response.
    """
    if (
        system.V_plus is None
        or system.V_minus is None
    ):
        raise ValueError(
            "RWA calculations require dipole_plus and "
            "dipole_minus in SpectroscopySystem"
        )

    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    return short_pulse_linear_rwa_response(
        L_super=system.L,
        V_plus=system.V_plus,
        V_minus=system.V_minus,
        observable=system.dipole,
        rho0=rho0,
        omega=omega,
        eta=eta,
        omega_L=pulse.omega_L,
        sigma=pulse.sigma,
        E0=pulse.E0,
        phase=pulse.phase,
        sign=sign,
    )


# ============================================================
# One-dimensional spectra
# ============================================================

def impulsive_linear_spectrum(
    system,
    omega,
    eta,
    rho0=None,
):
    """
    Full-interaction impulsive linear spectrum.
    """
    omega = np.asarray(
        omega,
        dtype=float,
    )

    if omega.ndim != 1:
        raise ValueError(
            "omega must be one-dimensional"
        )

    return np.array(
        [
            impulsive_linear_signal(
                system=system,
                omega=w,
                eta=eta,
                rho0=rho0,
            )
            for w in omega
        ],
        dtype=complex,
    )


def finite_pulse_linear_spectrum(
    system,
    pulse,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    Finite-Gaussian full-interaction linear spectrum for
    one selected field-frequency sign.
    """
    omega = np.asarray(
        omega,
        dtype=float,
    )

    if omega.ndim != 1:
        raise ValueError(
            "omega must be one-dimensional"
        )

    return np.array(
        [
            finite_pulse_linear_signal(
                system=system,
                pulse=pulse,
                omega=w,
                eta=eta,
                sign=sign,
                rho0=rho0,
            )
            for w in omega
        ],
        dtype=complex,
    )


def impulsive_linear_rwa_spectrum(
    system,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    Impulsive RWA linear spectrum.
    """
    omega = np.asarray(
        omega,
        dtype=float,
    )

    if omega.ndim != 1:
        raise ValueError(
            "omega must be one-dimensional"
        )

    return np.array(
        [
            impulsive_linear_rwa_signal(
                system=system,
                omega=w,
                eta=eta,
                sign=sign,
                rho0=rho0,
            )
            for w in omega
        ],
        dtype=complex,
    )


def short_pulse_linear_rwa_spectrum(
    system,
    pulse,
    omega,
    eta,
    sign=1,
    rho0=None,
):
    """
    Gaussian scalar-dressed RWA linear spectrum.
    """
    omega = np.asarray(
        omega,
        dtype=float,
    )

    if omega.ndim != 1:
        raise ValueError(
            "omega must be one-dimensional"
        )

    return np.array(
        [
            short_pulse_linear_rwa_signal(
                system=system,
                pulse=pulse,
                omega=w,
                eta=eta,
                sign=sign,
                rho0=rho0,
            )
            for w in omega
        ],
        dtype=complex,
    )
