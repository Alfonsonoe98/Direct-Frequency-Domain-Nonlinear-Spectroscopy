# PULSUS

**PULSUS**: **PU**lse-dressed **L**iouvillian **S**pectroscopy **U**sing **S**uperoperators.

PULSUS is a Python package for direct frequency-domain spectroscopy of open quantum systems using Liouville-space superoperators.

The current `v0.1` development version focuses on linear and third-order spectroscopy with time-independent Markovian Liouvillian dynamics, including finite Gaussian pulses, impulsive limits, rotating-wave approximations, arbitrary field-sign sectors, and complementary impulsive time-domain calculations.


## Core idea

For a field-free Liouvillian

\[
\mathcal{L},
\]

PULSUS evaluates frequency-domain propagation using resolvent actions

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


## Standard finite-pulse assumption

The standard finite-pulse formulation assumes that the state immediately before the pulse sequence is stationary under the field-free Liouvillian,

\[
\mathcal{L}|\rho_0\rangle\rangle = 0.
\]

For dissipative systems this is typically the stationary state of the Liouvillian.

PULSUS can also accept explicitly supplied initial density matrices for calculations where this assumption is not desired.


## Third-order spectroscopy

PULSUS currently provides four principal third-order models:

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


## Linear spectroscopy

The corresponding linear-response calculations are available through:

```python
pulsus.finite_pulse_linear_spectrum
pulsus.impulsive_linear_spectrum
pulsus.short_pulse_linear_rwa_spectrum
pulsus.impulsive_linear_rwa_spectrum
```


## Installation

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Run the test suite with:

```bash
pytest -q
```


## Basic usage

A spectroscopy model is defined using `SpectroscopySystem`:

```python
import pulsus

system = pulsus.SpectroscopySystem(
    H=H,
    dipole=mu,
    collapse_ops=collapse_ops,
)
```

The stationary state can be obtained from:

```python
rho_ss = system.stationary_state()
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

The returned complex array has shape:

```text
(len(omega3), len(omega1))
```

with convention

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


## Arbitrary field-sign sectors

PULSUS is not restricted to R and NR pathways.

Any third-order signature

\[
(s_1,s_2,s_3),
\qquad
s_j=\pm1,
\]

can be supplied explicitly.

For example:

```python
signature = (+1, +1, +1)

S = pulsus.finite_pulse_spectrum(
    system=system,
    pulse1=pulse,
    pulse2=pulse,
    pulse3=pulse,
    omega1=omega1,
    omega3=omega3,
    T=20.0,
    eta=0.02,
    signature=signature,
)
```

The aliases

```python
pathway="NR"
pathway="R"
```

remain available for the standard nonrephasing and rephasing sectors.


## Sign-reversed conjugation symmetry

For the tested coupled-dimer model, the field-sign sectors satisfy

\[
S_{(-s_1,-s_2,-s_3)}
(-\omega_1,-\omega_3)
=
S_{(s_1,s_2,s_3)}
(\omega_1,\omega_3)^*
\]

to floating-point precision.

This groups the eight possible signatures into four sign-reversed conjugate pairs.


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
    dipole_plus=mu_plus,
    dipole_minus=mu_minus,
)
```

PULSUS currently requires this decomposition to be supplied explicitly rather than inferring it automatically.


## Time-domain calculations

PULSUS also provides impulsive time-domain response functions.

For linear spectroscopy:

```python
R_t = pulsus.impulsive_linear_rwa_time_signal(
    system=system,
    times=times,
    sign=+1,
)
```

For third-order spectroscopy:

```python
R_t3_t1 = pulsus.impulsive_third_order_rwa_time_signal(
    system=system,
    t1=t1,
    t3=t3,
    T=20.0,
    pathway="NR",
)
```

These provide a complementary time-domain representation of the same impulsive Liouvillian dynamics used in the direct frequency-domain formulation.


## Pathway polarization in laboratory time

A selected RWA pathway can also be followed through the sequence of three impulsive interactions:

```python
P, components = pulsus.impulsive_pathway_polarization(
    system=system,
    times=times,
    pulse_times=(5.0, 15.0, 35.0),
    signature=(+1, -1, +1),
    return_components=True,
)
```

The resulting polarization is piecewise associated with the perturbative state active between pulses:

\[
P_{\mathrm{path}}(t)
=
\begin{cases}
0,
& t<\tau_1,\\
P^{(1)}(t),
& \tau_1\le t<\tau_2,\\
P^{(2)}(t),
& \tau_2\le t<\tau_3,\\
P^{(3)}(t),
& t\ge\tau_3.
\end{cases}
\]

This is the polarization of the **selected pathway**, not the total physical polarization from all perturbative orders and all field-sign sectors.


## Complex spectra and plotting

PULSUS retains the complete complex response.

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

Normalization is optional and disabled by default.


## Numerical implementation

Two-dimensional spectrum calculations exploit the separable dependence on the excitation and detection frequencies.

For a third-order response, the grid can be written schematically as

\[
S_{ji}=L_jR_i,
\]

where \(R_i\) contains the \(\omega_1\)-dependent operations and \(L_j\) contains the \(\omega_3\)-dependent operations.

This avoids repeating expensive matrix exponentials and resolvent solves independently at every pair

\[
(\omega_1,\omega_3).
\]

The optimized implementation reproduces the original point-by-point calculations to floating-point precision.


## Validation

PULSUS is tested against independently validated coupled-dimer calculations, frozen reference spectra, and explicit time-domain propagation.

For the finite-pulse development benchmark at

\[
\sigma=1.25,
\qquad
T=20,
\qquad
\eta=0.02,
\]

the production finite-pulse calculation gives complex spectral-shape discrepancies relative to an explicit finite-pulse time-domain ODE calculation of approximately

\[
0.61\%
\quad \text{(NR)}
\]

and

\[
0.64\%
\quad \text{(R)}.
\]

For comparison:

| Model | NR error | R error |
| --- | ---: | ---: |
| Finite-width PULSUS | 0.61% | 0.64% |
| Short-pulse RWA | 5.51% | 5.81% |
| Impulsive RWA | 9.60% | 9.83% |

These values characterize this development benchmark and are not universal accuracy bounds.


## Time-domain versus frequency-domain validation

For impulsive linear RWA response, numerical time propagation followed by an FFT agrees with the direct frequency-domain resolvent calculation with a relative complex error of approximately

\[
0.0035\%.
\]

For impulsive third-order NR response, two-dimensional time propagation followed by a 2D FFT agrees with the direct frequency-domain result with a relative complex error of approximately

\[
0.015\%.
\]

These comparisons independently validate the time/frequency sign and transform conventions used in PULSUS.


## Test suite

The current test suite contains **63 passing unit and regression tests** covering:

- Liouville-space vectorization and superoperators
- Lindblad Liouvillian construction
- stationary states
- resolvent actions
- Gaussian pulse dressing
- pulse-integral identities
- linear response
- third-order response
- R and NR sign conventions
- arbitrary field-sign signatures
- finite-pulse, impulsive, short-pulse, and RWA models
- two-dimensional spectrum construction
- optimized versus point-by-point evaluation
- plotting and array orientation
- frozen finite-pulse regression data
- impulsive time-domain propagation
- time-domain pathway polarization


## Examples

The `examples/` directory contains the current development and validation examples.

### Main examples

`01_finite_pulse_dimer.py`

Builds the coupled dimer and calculates finite-pulse NR and R spectra through the PULSUS API.

`03_reproduce_figure6.py`

Compares finite-pulse, short-pulse RWA, and impulsive RWA spectra.

`06_all_third_order_signatures.py`

Calculates all eight third-order field-sign sectors and checks the sign-reversed conjugation relation.

`09_linear_td_fft_vs_fd.py`

Validates impulsive linear time-domain propagation and FFT against the direct frequency-domain resolvent calculation.

`10_third_order_td_fft_vs_fd.py`

Validates third-order NR time-domain propagation and a 2D FFT against the direct frequency-domain calculation.

`11_impulsive_polarization_traces.py`

Visualizes the selected impulsive pathway in laboratory time and separates first-, second-, and third-order polarization stages.


### Validation and diagnostic examples

`02_compare_notebook04.py`

Compares the production finite-pulse implementation with frozen explicit-ODE validation data.

`04_linear_response_comparison.py`

Compares full-interaction and RWA linear spectra.

`05_linear_rwa_quadratures.py`

Examines real and imaginary quadratures of the linear full-interaction and RWA responses.

`07_full_vs_rwa_signatures.py`

Compares finite-pulse full-interaction and short-pulse RWA results across third-order signatures.

`08_approximation_ladder.py`

Separates the effects of the full-interaction, RWA, short-pulse, and impulsive approximations.


## Current scope

The current `v0.1` development version assumes:

- finite-dimensional quantum systems
- time-independent Liouvillian dynamics
- Markovian Lindblad evolution
- electric-dipole light-matter coupling
- Gaussian finite-duration pulses for the finite-pulse frequency-domain formulation
- temporally ordered, separated pulses
- linear and third-order response
- full-interaction and RWA calculations
- arbitrary third-order field-sign signatures
- impulsive time-domain propagation

Finite-duration driven time-domain propagation is not currently part of the public production API.


## Public API philosophy

The principal user-facing objects and functions are:

```python
pulsus.SpectroscopySystem
pulsus.GaussianPulse

pulsus.finite_pulse_spectrum
pulsus.impulsive_spectrum
pulsus.short_pulse_rwa_spectrum
pulsus.impulsive_rwa_spectrum

pulsus.finite_pulse_linear_spectrum
pulsus.impulsive_linear_spectrum
pulsus.short_pulse_linear_rwa_spectrum
pulsus.impulsive_linear_rwa_spectrum

pulsus.impulsive_linear_time_signal
pulsus.impulsive_third_order_rwa_time_signal
pulsus.impulsive_pathway_polarization

pulsus.plot_spectrum
```

Additional lower-level response, Liouville-space, resolvent, and pulse-dressing functions remain available for advanced calculations and validation.


## Repository structure

```text
.
├── src/
│   └── pulsus/          PULSUS source code
├── tests/               unit and regression tests
├── examples/            usage and validation examples
├── notebooks/           numerical validation notebooks
├── figures/             manuscript figures
├── results/             benchmark results
├── main.tex             manuscript / technical development
├── references.bib       bibliography
└── literature_notes.md  literature notes
```

PULSUS is under active development.
