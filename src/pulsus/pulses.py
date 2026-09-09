import numpy as np
from scipy.linalg import expm


def gaussian_prefactor(sigma, E0=1.0, phase=0.0, sign=1):
    """
    Complex prefactor for one signed-frequency component
    of a Gaussian pulse.

    C exp(i s phase),

    where

        C = (E0 / 2) sqrt(2 pi) sigma.

    Parameters
    ----------
    sigma : float
        Gaussian temporal width.
    E0 : float, optional
        Peak field amplitude.
    phase : float, optional
        Carrier phase.
    sign : {+1, -1}, optional
        Signed-frequency component.

    Returns
    -------
    complex
        Gaussian spectral prefactor.
    """
    if sign not in (-1, 1):
        raise ValueError("sign must be +1 or -1")

    C = 0.5 * E0 * np.sqrt(2.0 * np.pi) * sigma

    return C * np.exp(1j * sign * phase)


def _gaussian_dressing(
    L_super,
    detuning,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    Construct a Gaussian Liouvillian dressing operator

        F = C exp(i s phase)
            exp[(sigma^2 / 2) (L + i detuning I)^2].
    """
    L_super = np.asarray(L_super, dtype=complex)

    n = L_super.shape[0]
    I = np.eye(n, dtype=complex)

    A = L_super + 1j * detuning * I

    prefactor = gaussian_prefactor(
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )

    return prefactor * expm(
        0.5 * sigma**2 * (A @ A)
    )


def pulse_dressing_1(
    L_super,
    omega1,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    First-pulse Gaussian dressing operator.

        F1^(s)(omega1)
        = C exp(i s phase)
          exp[(sigma^2 / 2)
              (L + i(omega1 - s omega_L) I)^2].
    """
    detuning = omega1 - sign * omega_L

    return _gaussian_dressing(
        L_super,
        detuning=detuning,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )


def pulse_dressing_2(
    L_super,
    omega1,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    Second-pulse Gaussian dressing operator.

        F2^(s)(omega1)
        = C exp(i s phase)
          exp[(sigma^2 / 2)
              (L + i(omega1 + s omega_L) I)^2].
    """
    detuning = omega1 + sign * omega_L

    return _gaussian_dressing(
        L_super,
        detuning=detuning,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )


def pulse_dressing_3(
    L_super,
    omega3,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
):
    """
    Third-pulse Gaussian dressing operator.

        F3^(s)(omega3)
        = C exp(i s phase)
          exp[(sigma^2 / 2)
              (L + i(omega3 - s omega_L) I)^2].
    """
    detuning = omega3 - sign * omega_L

    return _gaussian_dressing(
        L_super,
        detuning=detuning,
        sigma=sigma,
        E0=E0,
        phase=phase,
        sign=sign,
    )
