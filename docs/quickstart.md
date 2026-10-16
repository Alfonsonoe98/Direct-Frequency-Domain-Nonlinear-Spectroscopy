# PULSUS quick start

## Install

PULSUS requires **Python 3.10 or newer**.

```bash
pip install pulsus-spec
```

The package name is `pulsus-spec`; import it in Python as `pulsus`. This guide describes version **0.2.0**. To install that version specifically, use `pip install "pulsus-spec==0.2.0"`.

## Your first spectrum

Copy this complete example into a notebook, or save it as `first_spectrum.py` and run `python first_spectrum.py`. It uses only the installed package and its dependencies. You do not need to clone the repository or use a particular folder. The program displays a figure without saving files.

```python
import matplotlib.pyplot as plt
import numpy as np
import pulsus

# Use consistent model units with hbar = 1.
# Build two coupled two-level sites.
identity = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)
a = np.kron(identity, lower)
b = np.kron(lower, identity)
H = (
    1.0 * (a.conj().T @ a)
    + 1.2 * (b.conj().T @ b)
    + 0.1 * (a.conj().T @ b + b.conj().T @ a)
)
dipole_plus = a.conj().T + b.conj().T
dipole_minus = dipole_plus.conj().T
dipole = dipole_plus + dipole_minus

gamma = 0.05
nbar = 0.05
collapse_ops = [
    np.sqrt(gamma * (nbar + 1)) * a,
    np.sqrt(gamma * nbar) * a.conj().T,
    np.sqrt(gamma * (nbar + 1)) * b,
    np.sqrt(gamma * nbar) * b.conj().T,
]
system = pulsus.SpectroscopySystem(
    H=H,
    dipole=dipole,
    collapse_ops=collapse_ops,
    hbar=1.0,
    dipole_plus=dipole_plus,
    dipole_minus=dipole_minus,
)
rho_ss = system.stationary_state()

# Define the pulse and frequency axes, then calculate the spectrum.
pulse = pulsus.GaussianPulse(omega_L=1.10, sigma=1.25, E0=1.0)
omega1 = np.linspace(0.7, 1.5, 31)
omega3 = np.linspace(0.7, 1.5, 31)
S = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    pathway="NR",
    rho0=rho_ss,
)

# Display the normalized imaginary component.
pulsus.plot_spectrum(
    omega1=omega1,
    omega3=omega3,
    spectrum=S,
    component="imag",
    normalize=True,
    title="Finite-pulse nonrephasing dimer spectrum",
    xlabel=r"$\omega_1$ (inverse time units)",
    ylabel=r"$\omega_3$ (inverse time units)",
)
plt.show()
```

The basic workflow is:

- `SpectroscopySystem` combines the Hamiltonian, dipole and collapse operators.
- `GaussianPulse` defines the carrier frequency, Gaussian field width and amplitude.
- `finite_pulse_spectrum` calculates the complex third-order response on two frequency axes.
- `plot_spectrum` displays a selected component of that response.

`stationary_state()` returns the initial density matrix; pass it explicitly as `rho0`. The raising and lowering dipole components allow this system to be reused with RWA models. The example uses consistent units with `hbar=1`: frequencies and rates are inverse time, while `sigma` and `T` are time. `sigma` is the standard deviation of the Gaussian **field** envelope.

`eta` is a positive convergence parameter for the causal frequency transform. It also adds broadening beyond physical dissipation. The example uses `eta=0.02` for illustration; for quantitative work, reduce `eta` and refine the frequency grid together until the features of interest are stable.

`S[j, i]` corresponds to `omega3[j]` and `omega1[i]`, so the result has shape `(31, 31)`. The plotting call displays the imaginary component divided by its maximum absolute value. It does not change the complex response stored in `S`.

## Choose a model

RWA means the **rotating-wave approximation**, which keeps resonant raising and lowering interactions. Five functions calculate third-order spectra:

| Function | Model |
|---|---|
| `finite_pulse_spectrum` (A) | Full interaction with finite Gaussian pulse dressing |
| `finite_pulse_rwa_spectrum` (B) | Finite Gaussian pulse dressing with RWA interactions |
| `short_pulse_rwa_spectrum` (C) | RWA molecular response multiplied by scalar Gaussian spectral envelopes |
| `impulsive_rwa_spectrum` (D) | RWA molecular response without pulse filtering or field prefactors |
| `impulsive_spectrum` | Full-interaction molecular response without pulses or an R/NR selector |

A–D form an approximation ladder. The full impulsive model is an additional reference.

To try B or C, change `pulsus.finite_pulse_spectrum` in the example to `pulsus.finite_pulse_rwa_spectrum` or `pulsus.short_pulse_rwa_spectrum` and keep every argument. For D, omit the pulse arguments:

```python
S = pulsus.impulsive_rwa_spectrum(
    system=system,
    omega1=omega1, omega3=omega3,
    T=20.0, eta=0.02,
    pathway="NR", rho0=rho_ss,
)
```

For the full impulsive model, use `pulsus.impulsive_spectrum` and also omit `pathway`. Impulsive responses have no Gaussian field factors, so their absolute amplitudes are not on the finite-pulse scale.

For linear spectroscopy, use `finite_pulse_linear_spectrum`, `short_pulse_linear_rwa_spectrum`, `impulsive_linear_spectrum` or `impulsive_linear_rwa_spectrum`. These return a complex array on one frequency axis.

Read the [function guide](api.md) to choose grid, single-coordinate and selected-coordinate calculations. The [user manual](user_manual/main.tex) explains system construction, units, pulses and model conventions.
