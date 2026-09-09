import numpy as np


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
