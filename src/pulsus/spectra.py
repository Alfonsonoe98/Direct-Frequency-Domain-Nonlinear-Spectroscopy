import numpy as np
from .response import (
    finite_pulse_signal,
    short_pulse_rwa_signal,
    impulsive_rwa_signal,
)
from .response import (
    finite_pulse_signal,
    impulsive_signal,
    short_pulse_rwa_signal,
    impulsive_rwa_signal,
)

def impulsive_spectrum(
    system,
    omega1,
    omega3,
    T,
    eta,
    rho0=None,
):
    """
    Evaluate the impulsive full-interaction third-order response
    on a two-dimensional frequency grid.
    """
    return response_grid(
        response_function=impulsive_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=T,
        eta=eta,
        rho0=rho0,
    )

def response_grid(
    response_function,
    omega1,
    omega3,
    **kwargs,
):
    """
    Evaluate a frequency-domain response on a 2D spectral grid.

    Parameters
    ----------
    response_function : callable
        Function that accepts omega1 and omega3 as keyword arguments
        and returns a complex response value.
    omega1 : array_like
        Excitation-frequency grid.
    omega3 : array_like
        Detection-frequency grid.
    **kwargs
        Additional keyword arguments passed to response_function.

    Returns
    -------
    numpy.ndarray
        Complex response array with shape
        (len(omega3), len(omega1)).
    """
    omega1 = np.asarray(omega1, dtype=float)
    omega3 = np.asarray(omega3, dtype=float)

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError("omega1 and omega3 must be one-dimensional")

    spectrum = np.empty(
        (omega3.size, omega1.size),
        dtype=complex,
    )

    for j, w3 in enumerate(omega3):
        for i, w1 in enumerate(omega1):
            spectrum[j, i] = response_function(
                omega1=w1,
                omega3=w3,
                **kwargs,
            )

    return spectrum

def finite_pulse_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate a finite-pulse third-order spectrum on a 2D
    frequency grid.

    Parameters
    ----------
    system : SpectroscopySystem
        Open quantum system.
    pulse1, pulse2, pulse3 : GaussianPulse
        Three Gaussian laser pulses.
    omega1 : array_like
        Signed excitation-frequency grid.
    omega3 : array_like
        Detection-frequency grid.
    T : float
        Waiting time.
    eta : float
        Resolvent broadening parameter.
    pathway : {"NR", "R"}, optional
        Nonrephasing or rephasing field-sign sector.
    rho0 : array_like, optional
        Initial density matrix. If omitted, system.rho0 is used.

    Returns
    -------
    numpy.ndarray
        Complex third-order spectrum with shape

            (len(omega3), len(omega1)).
    """
    return response_grid(
        response_function=finite_pulse_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        T=T,
        eta=eta,
        pathway=pathway,
        rho0=rho0,
    )

def short_pulse_rwa_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate the short-pulse RWA third-order spectrum on a
    two-dimensional frequency grid.
    """
    return response_grid(
        response_function=short_pulse_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        pulse1=pulse1,
        pulse2=pulse2,
        pulse3=pulse3,
        T=T,
        eta=eta,
        pathway=pathway,
        rho0=rho0,
    )


def impulsive_rwa_spectrum(
    system,
    omega1,
    omega3,
    T,
    eta,
    pathway="NR",
    rho0=None,
):
    """
    Evaluate the impulsive RWA third-order spectrum on a
    two-dimensional frequency grid.
    """
    return response_grid(
        response_function=impulsive_rwa_signal,
        omega1=omega1,
        omega3=omega3,
        system=system,
        T=T,
        eta=eta,
        pathway=pathway,
        rho0=rho0,
    )
