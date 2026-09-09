import numpy as np

import pulsus
from dimer_model import build_dimer


def test_vec_unvec_roundtrip():
    rho = np.array(
        [[1.0, 2.0j],
         [-2.0j, 3.0]],
        dtype=complex,
    )

    rho_back = pulsus.unvec(pulsus.vec(rho), 2)

    assert np.allclose(rho_back, rho)


def test_hamiltonian_liouvillian_action():
    H = np.array(
        [[0.0, 0.2],
         [0.2, 1.0]],
        dtype=complex,
    )

    rho = np.array(
        [[0.7, 0.1j],
         [-0.1j, 0.3]],
        dtype=complex,
    )

    L_H = pulsus.hamiltonian_liouvillian(H)

    result = pulsus.unvec(
        L_H @ pulsus.vec(rho),
        2,
    )

    reference = -1j * (H @ rho - rho @ H)

    assert np.allclose(result, reference)


def test_interaction_superoperator_action():
    mu = np.array(
        [[0.0, 1.0],
         [1.0, 0.0]],
        dtype=complex,
    )

    rho = np.array(
        [[1.0, 0.0],
         [0.0, 0.0]],
        dtype=complex,
    )

    V = pulsus.interaction_superoperator(mu)

    result = pulsus.unvec(
        V @ pulsus.vec(rho),
        2,
    )

    reference = 1j * (mu @ rho - rho @ mu)

    assert np.allclose(result, reference)


def test_full_liouvillian_matches_dimer_reference():
    model = build_dimer()

    collapse_ops = [
        model["L_a_down"],
        model["L_a_up"],
        model["L_b_down"],
        model["L_b_up"],
    ]

    L = pulsus.liouvillian(
        model["H_S"],
        collapse_ops=collapse_ops,
        hbar=model["hbar"],
    )

    assert np.allclose(L, model["L_super"])


def test_stationary_state():
    model = build_dimer()

    rho_ss = pulsus.stationary_state(
        model["L_super"],
        model["d"],
    )

    assert np.allclose(np.trace(rho_ss), 1.0)

    residual = model["L_super"] @ pulsus.vec(rho_ss)

    assert np.linalg.norm(residual) < 1e-10
