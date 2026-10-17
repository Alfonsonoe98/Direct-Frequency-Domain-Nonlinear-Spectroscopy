# PULSUS

**PU**lse-dressed **L**iouvillian **S**pectroscopy **U**sing **S**uperoperators.

PULSUS calculates linear and third-order spectra of open quantum systems directly
in the frequency domain. It supports finite Gaussian pulses, full and rotating-wave
interactions, impulsive limits, and selected spectral coordinates.

## Install

Requires **Python 3.10 or newer**:

```bash
pip install pulsus-spec
```

Import the package as `pulsus`. The first example below works with **0.2.0 or
newer**. The linear plotter and optional figure export require **0.3.0 or
newer**. To install from a cloned checkout instead of PyPI, run this from the
repository root:

```bash
python -m pip install -e .
```

## Documentation

- [Quick start](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/quickstart.md): install, run a first calculation, and choose a model.
- [User manual (PDF)](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/user_manual/PULSUS_User_Manual.pdf): systems, pulses, units, and complete tutorials.
- [API reference](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/api.md): public functions and their arguments.
- [Method and validation notes](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/method_notes.md): scientific conventions, implementation details, and benchmark results.
- [Examples](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/tree/main/examples/manual) and [documentation source/build instructions](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/README.md). Overleaf is optional.

## Your first spectrum

Copy this complete example into a notebook or Python file. It uses only the
installed package and its dependencies, displays a finite-pulse nonrephasing
spectrum of two coupled sites, and saves no files.

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

![Finite-pulse nonrephasing dimer spectrum](https://raw.githubusercontent.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/main/docs/user_manual/figures/first_nr_spectrum.png)

`S[j, i]` corresponds to `omega3[j]` and `omega1[i]`. All quantities use
consistent model units with `hbar=1`; `sigma` is the Gaussian **field-envelope**
standard deviation. `eta` ensures convergence of the causal transform and adds
broadening: reduce it and refine the frequency grid until the result is stable.

## Available models

The rotating-wave approximation (RWA) and impulsive approximation can be chosen
separately. The principal third-order grid functions are:

| Function | Pulse treatment | Interaction |
|---|---|---|
| `finite_pulse_spectrum` | Finite Gaussian pulses | Full |
| `finite_pulse_rwa_spectrum` | Finite Gaussian pulses | RWA |
| `short_pulse_rwa_spectrum` | Gaussian spectral envelopes | RWA |
| `impulsive_rwa_spectrum` | Impulsive | RWA |
| `impulsive_spectrum` | Impulsive | Full |

For linear spectra, use `finite_pulse_linear_spectrum`,
`short_pulse_linear_rwa_spectrum`, `impulsive_linear_rwa_spectrum`, or
`impulsive_linear_spectrum`. Functions return complex numerical responses;
plotting is a separate step. See the [model guide](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/quickstart.md#choose-a-model)
for the arguments to keep or omit when switching models.

## Plotting and optional export in 0.3.0

Using `system`, `pulse`, and `rho_ss` from the first example:

```python
omega = np.linspace(0.7, 1.5, 401)
S_linear = pulsus.finite_pulse_linear_spectrum(
    system=system, pulse=pulse, omega=omega, eta=0.02, rho0=rho_ss,
)
fig, ax, line = pulsus.plot_linear_spectrum(
    omega, S_linear, component="imag", normalize=True,
    title="Finite-pulse linear spectrum",
)
plt.show()
```

`plot_linear_spectrum` draws a line without markers by default. Both plotters
support `component="real"`, `"imag"`, or `"abs"`, optional normalization,
`savepath="spectrum.png"` or `"spectrum.pdf"`, and `dpi=300`. They return
editable Matplotlib objects without calling `plt.show()` or closing the figure.
Omitting `savepath` saves nothing. You can modify `ax` and then call
`fig.savefig(...)` and `plt.show()`.
See the [plotting reference](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/api.md#5-plot-linear-and-third-order-spectra)
for overlays and all options.

## Physical scope

The formulation uses finite-dimensional systems, electric-dipole coupling, and a
time-independent Markovian Lindblad Liouvillian. Standard finite-pulse
calculations assume a stationary initial state and temporally ordered,
separated Gaussian pulses. RWA calculations require explicit raising and lowering
dipole components. Impulsive time-domain response functions are also available;
finite-duration driven time propagation is not part of the public API.

## Development

From a cloned repository:

```bash
python -m pip install -e ".[dev]"
pytest -q
```

The suite includes 108 unit and regression tests. See the
[method and validation notes](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/method_notes.md#validation)
for scientific benchmarks and the [changelog](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/CHANGELOG.md)
for release changes. Maintainers can follow the
[release checklist](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/docs/release_checklist.md)
to publish a new version.

## Citation and license

The archived 0.2.0 release can be cited as:

> Castillo-Gonzalez, Alfonso. *PULSUS: Pulse-dressed Liouvillian Spectroscopy Using Superoperators*. Version 0.2.0, Zenodo, 2026. [doi:10.5281/zenodo.23253040](https://doi.org/10.5281/zenodo.23253040).

See [CITATION.cff](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/CITATION.cff) for citation metadata.
PULSUS is distributed under the [BSD 3-Clause license](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/blob/main/LICENSE).
