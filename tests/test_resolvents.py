import numpy as np

import pulsus
from dimer_model import build_dimer


def test_resolvent_matrix_defining_equation():
    model = build_dimer()

    L = model["L_super"]
    omega = 1.0
    eta = 0.02

    G = pulsus.resolvent_matrix(
        L,
        omega=omega,
        eta=eta,
    )

    n = L.shape[0]
    I = np.eye(n, dtype=complex)

    A = (eta - 1j * omega) * I - L

    assert np.allclose(A @ G, I)


def test_resolvent_action():
    model = build_dimer()

    L = model["L_super"]
    rhs = pulsus.vec(model["rho0"])

    omega = 1.0
    eta = 0.02

    x = pulsus.resolvent_action(
        L,
        omega=omega,
        rhs=rhs,
        eta=eta,
    )

    n = L.shape[0]
    I = np.eye(n, dtype=complex)

    A = (eta - 1j * omega) * I - L

    assert np.allclose(A @ x, rhs)


def test_resolvent_action_matches_matrix():
    model = build_dimer()

    L = model["L_super"]
    rhs = pulsus.vec(model["rho0"])

    omega = 1.0
    eta = 0.02

    G = pulsus.resolvent_matrix(L, omega, eta)

    x_matrix = G @ rhs
    x_action = pulsus.resolvent_action(
        L,
        omega,
        rhs,
        eta,
    )

    assert np.allclose(x_matrix, x_action)
