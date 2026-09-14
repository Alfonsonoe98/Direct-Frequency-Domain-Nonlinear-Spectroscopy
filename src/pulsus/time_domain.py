import numpy as np
from scipy.linalg import expm

from .liouville import vec
from .response import observable_bra, resolve_signature


def _validate_time_axis(times, name):
    times = np.asarray(
        times,
        dtype=float,
    )

    if times.ndim != 1:
        raise ValueError(
            f"{name} must be one-dimensional"
        )

    if np.any(times < 0.0):
        raise ValueError(
            f"{name} must contain non-negative times"
        )

    return times


def _initial_state(system, rho0):
    if rho0 is None:
        rho0 = system.rho0

    if rho0 is None:
        raise ValueError(
            "An initial state must be supplied either through "
            "system.rho0 or the rho0 argument"
        )

    return vec(rho0)


def _require_rwa(system):
    if system.V_plus is None or system.V_minus is None:
        raise ValueError(
            "RWA calculations require dipole_plus and "
            "dipole_minus in SpectroscopySystem"
        )


def impulsive_linear_time_signal(
    system,
    times,
    rho0=None,
):
    """
    Evaluate the impulsive full-V linear response.

    R^(1)(t) = <mu| exp(L t) V |rho0>
    """
    times = _validate_time_axis(
        times,
        "times",
    )

    state0 = _initial_state(
        system,
        rho0,
    )

    state = (
        system.V
        @ state0
    )

    bra = observable_bra(
        system.dipole
    )

    signal = np.empty(
        len(times),
        dtype=complex,
    )

    for i, t in enumerate(times):
        signal[i] = (
            bra
            @ expm(system.L * t)
            @ state
        )

    return signal


def impulsive_linear_rwa_time_signal(
    system,
    times,
    sign=1,
    rho0=None,
):
    """
    Evaluate one impulsive linear RWA field-sign sector.
    """
    _require_rwa(system)

    if sign not in (-1, 1):
        raise ValueError(
            "sign must be +1 or -1"
        )

    times = _validate_time_axis(
        times,
        "times",
    )

    state0 = _initial_state(
        system,
        rho0,
    )

    V = (
        system.V_plus
        if sign == 1
        else system.V_minus
    )

    state = (
        V
        @ state0
    )

    bra = observable_bra(
        system.dipole
    )

    signal = np.empty(
        len(times),
        dtype=complex,
    )

    for i, t in enumerate(times):
        signal[i] = (
            bra
            @ expm(system.L * t)
            @ state
        )

    return signal


def _third_order_time_grid(
    system,
    t1,
    t3,
    T,
    V1,
    V2,
    V3,
    rho0,
):
    t1 = _validate_time_axis(
        t1,
        "t1",
    )

    t3 = _validate_time_axis(
        t3,
        "t3",
    )

    if T < 0.0:
        raise ValueError(
            "T must be non-negative"
        )

    state0 = _initial_state(
        system,
        rho0,
    )

    bra = observable_bra(
        system.dipole
    )

    first = (
        V1
        @ state0
    )

    right = np.column_stack(
        [
            expm(system.L * time)
            @ first
            for time in t1
        ]
    )

    right = (
        expm(system.L * T)
        @ (
            V2
            @ right
        )
    )

    right = (
        V3
        @ right
    )

    signal = np.empty(
        (
            len(t3),
            len(t1),
        ),
        dtype=complex,
    )

    for j, time in enumerate(t3):
        signal[j, :] = (
            bra
            @ expm(system.L * time)
            @ right
        )

    return signal


def impulsive_third_order_time_signal(
    system,
    t1,
    t3,
    T,
    rho0=None,
):
    """
    Evaluate the full-V impulsive third-order response.

    R^(3)(t3,T,t1)
      = <mu|
        exp(L t3) V
        exp(L T) V
        exp(L t1) V
        |rho0>.
    """
    return _third_order_time_grid(
        system=system,
        t1=t1,
        t3=t3,
        T=T,
        V1=system.V,
        V2=system.V,
        V3=system.V,
        rho0=rho0,
    )


def impulsive_third_order_rwa_time_signal(
    system,
    t1,
    t3,
    T,
    pathway="NR",
    signature=None,
    rho0=None,
):
    """
    Evaluate one impulsive third-order RWA field-sign sector.
    """
    _require_rwa(system)

    s1, s2, s3 = resolve_signature(
        pathway=pathway,
        signature=signature,
    )

    V1 = (
        system.V_plus
        if s1 == 1
        else system.V_minus
    )

    V2 = (
        system.V_plus
        if s2 == 1
        else system.V_minus
    )

    V3 = (
        system.V_plus
        if s3 == 1
        else system.V_minus
    )

    return _third_order_time_grid(
        system=system,
        t1=t1,
        t3=t3,
        T=T,
        V1=V1,
        V2=V2,
        V3=V3,
        rho0=rho0,
    )


def impulsive_linear_polarization(
    system,
    times,
    pulse_time=0.0,
    sign=None,
    rho0=None,
):
    """
    Evaluate an impulsive linear polarization trace P^(1)(t).

    The response is zero before pulse_time.

    If sign is None, the complete interaction V is used.
    If sign is +1 or -1, the corresponding RWA interaction
    V_plus or V_minus is used.
    """
    times = _validate_time_axis(
        times,
        "times",
    )

    if pulse_time < 0.0:
        raise ValueError(
            "pulse_time must be non-negative"
        )

    if sign is not None:
        if sign not in (-1, 1):
            raise ValueError(
                "sign must be +1, -1, or None"
            )

        _require_rwa(system)

    state0 = _initial_state(
        system,
        rho0,
    )

    if sign is None:
        V = system.V
    elif sign == 1:
        V = system.V_plus
    else:
        V = system.V_minus

    state = (
        V
        @ state0
    )

    bra = observable_bra(
        system.dipole
    )

    polarization = np.zeros(
        len(times),
        dtype=complex,
    )

    mask = (
        times >= pulse_time
    )

    for i in np.where(mask)[0]:
        delay = (
            times[i]
            - pulse_time
        )

        polarization[i] = (
            bra
            @ expm(system.L * delay)
            @ state
        )

    return polarization


def impulsive_third_order_polarization(
    system,
    times,
    pulse_times,
    pathway=None,
    signature=None,
    rho0=None,
):
    """
    Evaluate a third-order impulsive polarization trace P^(3)(t)
    for one ordered three-pulse sequence.

    pulse_times = (tau1, tau2, tau3)

    The third-order contribution is zero before tau3.

    If pathway and signature are both None, the complete
    interaction V is used at all three interactions.

    If pathway or signature is supplied, the corresponding RWA
    field-sign sector is evaluated. Examples:

        pathway="NR"
        pathway="R"
        signature=(+1, -1, +1)
    """
    times = _validate_time_axis(
        times,
        "times",
    )

    pulse_times = np.asarray(
        pulse_times,
        dtype=float,
    )

    if pulse_times.shape != (3,):
        raise ValueError(
            "pulse_times must contain exactly "
            "(tau1, tau2, tau3)"
        )

    tau1, tau2, tau3 = pulse_times

    if np.any(pulse_times < 0.0):
        raise ValueError(
            "pulse_times must be non-negative"
        )

    if not (
        tau1 <= tau2 <= tau3
    ):
        raise ValueError(
            "pulse_times must satisfy "
            "tau1 <= tau2 <= tau3"
        )

    state0 = _initial_state(
        system,
        rho0,
    )

    use_rwa = (
        pathway is not None
        or signature is not None
    )

    if use_rwa:
        _require_rwa(system)

        if pathway is None:
            pathway = "NR"

        s1, s2, s3 = resolve_signature(
            pathway=pathway,
            signature=signature,
        )

        V1 = (
            system.V_plus
            if s1 == 1
            else system.V_minus
        )

        V2 = (
            system.V_plus
            if s2 == 1
            else system.V_minus
        )

        V3 = (
            system.V_plus
            if s3 == 1
            else system.V_minus
        )

    else:
        V1 = system.V
        V2 = system.V
        V3 = system.V

    t1 = (
        tau2
        - tau1
    )

    T = (
        tau3
        - tau2
    )

    state = (
        V1
        @ state0
    )

    state = (
        expm(system.L * t1)
        @ state
    )

    state = (
        V2
        @ state
    )

    state = (
        expm(system.L * T)
        @ state
    )

    state = (
        V3
        @ state
    )

    bra = observable_bra(
        system.dipole
    )

    polarization = np.zeros(
        len(times),
        dtype=complex,
    )

    mask = (
        times >= tau3
    )

    for i in np.where(mask)[0]:
        t3 = (
            times[i]
            - tau3
        )

        polarization[i] = (
            bra
            @ expm(system.L * t3)
            @ state
        )

    return polarization

def impulsive_pathway_polarization(
    system,
    times,
    pulse_times,
    pathway="NR",
    signature=None,
    rho0=None,
    return_components=False,
):
    """
    Evaluate the piecewise polarization history of one ordered
    impulsive RWA pathway.

    The pathway evolves through three stages:

        tau1 <= t < tau2:
            first-order coherence

        tau2 <= t < tau3:
            second-order state

        t >= tau3:
            third-order polarization

    Parameters
    ----------
    system
        SpectroscopySystem with V_plus and V_minus.

    times : array_like
        Laboratory-time axis.

    pulse_times : tuple
        (tau1, tau2, tau3).

    pathway : {"NR", "R"}
        Convenience field-sign alias.

    signature : tuple, optional
        Explicit (s1, s2, s3) field-sign signature.

    rho0 : ndarray, optional
        Initial density matrix.

    return_components : bool
        If True, also return the separate first-, second-,
        and third-stage traces.

    Returns
    -------
    polarization : ndarray
        Piecewise complex polarization history.

    components : dict, optional
        Returned only when return_components=True.
    """
    _require_rwa(system)

    times = _validate_time_axis(
        times,
        "times",
    )

    pulse_times = np.asarray(
        pulse_times,
        dtype=float,
    )

    if pulse_times.shape != (3,):
        raise ValueError(
            "pulse_times must contain exactly "
            "(tau1, tau2, tau3)"
        )

    tau1, tau2, tau3 = pulse_times

    if np.any(pulse_times < 0.0):
        raise ValueError(
            "pulse_times must be non-negative"
        )

    if not (
        tau1 <= tau2 <= tau3
    ):
        raise ValueError(
            "pulse_times must satisfy "
            "tau1 <= tau2 <= tau3"
        )

    s1, s2, s3 = resolve_signature(
        pathway=pathway,
        signature=signature,
    )

    V1 = (
        system.V_plus
        if s1 == 1
        else system.V_minus
    )

    V2 = (
        system.V_plus
        if s2 == 1
        else system.V_minus
    )

    V3 = (
        system.V_plus
        if s3 == 1
        else system.V_minus
    )

    state0 = _initial_state(
        system,
        rho0,
    )

    bra = observable_bra(
        system.dipole
    )

    # ========================================================
    # State immediately after pulse 1
    # ========================================================

    state1 = (
        V1
        @ state0
    )

    # ========================================================
    # State immediately after pulse 2
    # ========================================================

    state1_at_tau2 = (
        expm(
            system.L
            * (tau2 - tau1)
        )
        @ state1
    )

    state2 = (
        V2
        @ state1_at_tau2
    )

    # ========================================================
    # State immediately after pulse 3
    # ========================================================

    state2_at_tau3 = (
        expm(
            system.L
            * (tau3 - tau2)
        )
        @ state2
    )

    state3 = (
        V3
        @ state2_at_tau3
    )

    # ========================================================
    # Piecewise laboratory-time traces
    # ========================================================

    P1 = np.zeros(
        len(times),
        dtype=complex,
    )

    P2 = np.zeros(
        len(times),
        dtype=complex,
    )

    P3 = np.zeros(
        len(times),
        dtype=complex,
    )

    mask1 = (
        (times >= tau1)
        & (times < tau2)
    )

    mask2 = (
        (times >= tau2)
        & (times < tau3)
    )

    mask3 = (
        times >= tau3
    )

    for i in np.where(mask1)[0]:
        delay = (
            times[i]
            - tau1
        )

        P1[i] = (
            bra
            @ expm(
                system.L * delay
            )
            @ state1
        )

    for i in np.where(mask2)[0]:
        delay = (
            times[i]
            - tau2
        )

        P2[i] = (
            bra
            @ expm(
                system.L * delay
            )
            @ state2
        )

    for i in np.where(mask3)[0]:
        delay = (
            times[i]
            - tau3
        )

        P3[i] = (
            bra
            @ expm(
                system.L * delay
            )
            @ state3
        )

    polarization = (
        P1
        + P2
        + P3
    )

    if return_components:
        return polarization, {
            "first_order": P1,
            "second_order": P2,
            "third_order": P3,
        }

    return polarization
