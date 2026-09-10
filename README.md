# PULSUS

**PULSUS**: **PU**lse-dressed **L**iouvillian **S**pectroscopy **U**sing **S**uperoperators.

PULSUS is a Python package for direct frequency-domain nonlinear spectroscopy of open quantum systems using Liouville-space superoperators.

The current development version focuses on third-order spectroscopy with time-independent Markovian Liouvillian dynamics and Gaussian laser pulses.

## Core idea

For a field-free Liouvillian

\[
\mathcal{L},
\]

PULSUS evaluates frequency-domain propagation through resolvent actions

\[
G_\eta(\omega)
=
\left[
(\eta-i\omega)I-\mathcal{L}
\right]^{-1},
\]

without explicitly constructing the matrix inverse.

Finite Gaussian pulses are represented by Liouvillian dressing operators. For example,

\[
F_1^{(s)}
=
C e^{is\Phi}
\exp\left[
\frac{\sigma^2}{2}
\left(
\mathcal{L}
+i(\omega_1-s\omega_L)I
\right)^2
\right].
\]

The production finite-pulse implementation uses the physical field-free Liouvillian inside the pulse dressing,

\[
F_j = F_j(\mathcal{L}),
\]

while the resolvents may use a small positive broadening parameter \(\eta\).

## Implemented third-order models

PULSUS currently provides four related response models:

| Function | Pulse treatment | Light-matter interaction |
| --- | --- | --- |
| `finite_pulse_spectrum` | finite Gaussian pulses | full interaction |
| `impulsive_spectrum` | impulsive | full interaction |
| `short_pulse_rwa_spectrum` | Gaussian spectral envelopes | RWA |
| `impulsive_rwa_spectrum` | impulsive | RWA |

The impulsive approximation and the rotating-wave approximation are treated as distinct approximations.

In particular,

\[
\sigma \rightarrow 0
\]

for the full finite-pulse theory approaches the impulsive **full-interaction** response, not automatically the impulsive RWA response.

## Installation

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Run the test suite with:

```bash
pytest -v
```

## Basic usage

A spectroscopy model is defined using `SpectroscopySystem`:

```python
import pulsus

system = pulsus.SpectroscopySystem(
    H=H,
    dipole=mu,
    collapse_ops=collapse_ops,
    rho0=rho0,
)
```

A Gaussian pulse is specified by its carrier frequency, temporal width, amplitude, and phase:

```python
pulse = pulsus.GaussianPulse(
    omega_L=1.10,
    sigma=1.25,
    E0=1.0,
    phase=0.0,
)
```

A finite-pulse nonrephasing spectrum can then be evaluated directly:

```python
import numpy as np

omega1 = np.linspace(0.7, 1.5, 81)
omega3 = np.linspace(0.7, 1.5, 81)

S_NR = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    pathway="NR",
)
```

The returned complex array has shape

```text
(len(omega3), len(omega1))
```

with the convention

\[
S[j,i]
=
S(\omega_{3,j},\omega_{1,i}).
\]

## Rephasing and nonrephasing convention

PULSUS uses signed first-coherence frequencies.

For the nonrephasing sector,

\[
\omega_1>0,
\]

with field-sign signature

\[
(+,-,+).
\]

For the rephasing sector,

\[
\omega_1<0,
\]

with field-sign signature

\[
(-,+,+).
\]

For example:

```python
omega1_NR = np.linspace(0.7, 1.5, 81)
omega1_R = -omega1_NR
```

The rephasing spectrum may be displayed against the positive spectroscopic excitation coordinate `-omega1_R`, but the calculation itself uses signed negative first-coherence frequencies.

## Full interaction and RWA

Full-interaction calculations require only the total dipole operator

\[
\mu.
\]

RWA calculations additionally require a decomposition

\[
\mu = \mu_+ + \mu_-.
\]

For example:

```python
system = pulsus.SpectroscopySystem(
    H=H,
    dipole=mu,
    collapse_ops=collapse_ops,
    rho0=rho0,
    dipole_plus=mu_plus,
    dipole_minus=mu_minus,
)
```

PULSUS currently requires this decomposition to be supplied explicitly rather than inferring it automatically.

## Complex spectra and plotting

PULSUS retains the complete complex nonlinear response.

A spectrum can be plotted with:

```python
pulsus.plot_spectrum(
    omega1=omega1,
    omega3=omega3,
    spectrum=S_NR,
    component="imag",
    normalize=True,
)
```

Supported displayed components are:

- `"real"`
- `"imag"`
- `"abs"`

Normalization is optional and is disabled by default.

## Numerical implementation

Spectrum calculations exploit the separable dependence on the excitation and detection frequencies.

For a third-order response, the grid can be written schematically as

\[
S_{ji}=L_jR_i,
\]

where \(R_i\) contains the \(\omega_1\)-dependent operations and \(L_j\) contains the \(\omega_3\)-dependent operations.

This avoids repeating expensive matrix exponentials and resolvent solves independently at every pair \((\omega_1,\omega_3)\).

The optimized implementation preserves the results of the original point-by-point calculation to floating-point precision.

## Validation

The implementation is tested against independently validated coupled-dimer calculations and frozen reference spectra generated by the validation notebooks.

For the finite-pulse development benchmark at

\[
\sigma=1.25,\qquad
T=20,\qquad
\eta=0.02,
\]

the production PULSUS finite-pulse calculation gives complex spectral-shape discrepancies relative to an explicit finite-pulse time-domain ODE calculation of approximately

\[
0.61\% \quad \text{(NR)}
\]

and

\[
0.64\% \quad \text{(R)}.
\]

For comparison, under the same benchmark conditions:

| Model | NR error | R error |
| --- | ---: | ---: |
| Finite-width PULSUS | 0.61% | 0.64% |
| Short-pulse RWA | 5.51% | 5.81% |
| Impulsive RWA | 9.60% | 9.83% |

These values characterize this development benchmark and are not universal accuracy bounds.

The current test suite contains 38 passing unit and regression tests covering:

- Liouville-space vectorization and superoperators
- Lindblad Liouvillian construction
- stationary states
- resolvent actions
- Gaussian pulse dressing
- defining pulse-integral identities
- third-order response construction
- rephasing and nonrephasing sign conventions
- finite-pulse, impulsive, and RWA responses
- two-dimensional spectrum construction
- optimized versus point-by-point spectrum evaluation
- plotting and array orientation
- regression against frozen validation spectra

## Examples

The `examples/` directory currently contains:

- `01_finite_pulse_dimer.py`  
  Builds the coupled dimer and calculates finite-pulse NR and R spectra through the public PULSUS API.

- `02_compare_notebook04.py`  
  Compares the production finite-pulse implementation with frozen explicit-ODE validation data.

- `03_reproduce_figure6.py`  
  Reproduces the main finite-pulse, short-pulse RWA, and impulsive RWA comparison through the public PULSUS API.

## Current scope

The current development version assumes:

- finite-dimensional quantum systems
- time-independent Liouvillian dynamics
- Markovian Lindblad evolution
- electric-dipole light-matter coupling
- Gaussian finite-duration pulses
- temporally ordered, separated pulses
- third-order nonlinear response
- rephasing and nonrephasing field-sign sectors

PULSUS is under active development.

## Repository structure

```text
.
├── src/
│   └── pulsus/          PULSUS source code
├── tests/               unit and regression tests
├── examples/            public API examples
├── notebooks/           numerical validation notebooks
├── figures/             manuscript figures
├── results/             benchmark results
├── main.tex             manuscript / technical development
├── references.bib       bibliography
└── literature_notes.md  literature notes
```
