# PULSUS 0.2.0 public API reference

This reference describes the source in the repository's PULSUS 0.2.0 checkout. Install the distribution with `python -m pip install pulsus-spec` and import it with `import pulsus`. For a version-matched installation use `python -m pip install "pulsus-spec==0.2.0"`. Follow [the quick start](quickstart.md) for a complete first calculation.

All names below are available as `pulsus.NAME`. The root module imports **57 names**, while its `__all__` contains **55**. `finite_pulse_rwa_response` and `finite_pulse_rwa_signal` are importable explicitly but omitted from `__all__` in this source version. Thus `import pulsus` or explicit imports expose them; `from pulsus import *` omits them. This reference records the discrepancy without changing scientific code.

## Most commonly used functions

API means **application programming interface**: the public classes and functions you call from Python. Start with these:

| Step | Class or function | Use |
|---|---|---|
| 1 | `SpectroscopySystem` | Define the Hamiltonian, dipole, dissipation and optional state |
| 2 | `GaussianPulse` | Define carrier frequency, temporal width, amplitude and phase |
| 3 | `system.stationary_state()` | Obtain a stationary initial density matrix; pass it as `rho0` |
| 4 | `*_spectrum` | Calculate a linear frequency array or a 2D third-order grid |
| 5 | `plot_spectrum` | Display a 2D grid; use Matplotlib directly for linear spectra |
| Optional | `finite_pulse_points` | Evaluate selected frequency pairs without a Cartesian grid |

## Signal, spectrum, response, and time_signal

Choose a model prefix, such as `finite_pulse` or `impulsive_rwa`; the suffix specifies the interface:

| Suffix or function | Inputs | Result | Typical use |
|---|---|---|---|
| `*_spectrum` (linear) | system object and a 1D frequency axis | `(N,)` complex array | Plot a linear spectrum |
| `*_spectrum` (third order) | system object and two 1D frequency axes | `(N3, N1)` complex array | Plot a 2D spectrum |
| `*_signal` | system object and scalar frequency coordinate(s) | complex scalar | Evaluate one frequency or one 2D pair |
| `finite_pulse_points` | system object and equal-length paired frequency arrays | `(N,)` complex array | Evaluate a cut or selected coordinates |
| `*_response` | Liouvillian, interaction and detection matrices, plus scalar frequency coordinate(s) | complex scalar | Work directly with superoperators |
| `*_time_signal` | system object and delay-time arrays | complex vector or grid | Inspect impulsive evolution on time axes |
| `*_polarization` | system object and laboratory-time array | complex vector, optionally stage components | Place an impulsive trace relative to pulse arrival times |

For example, `finite_pulse_signal` returns one value at scalar `omega1` and `omega3`, while `finite_pulse_spectrum` returns all Cartesian combinations of two frequency arrays. These names do not add a detector model. Constructing a phased absorptive combination, heterodyne signal or intensity requires the intended physical convention; plotting `imag` or `abs` does not do that automatically.

## Conventions shared by the functions

- Use consistent units. `H/hbar`, `omega`, `omega1`, `omega3`, all carrier frequencies, decay rates, and `eta` have inverse-time units. `sigma`, `times`, `t1`, `t3`, `T`, and pulse arrival times have time units. Phases are radians. `hbar=1` is a consistent model-unit choice, not a unit conversion.
- Hilbert-space operators and density matrices are `(d, d)` arrays. Superoperators and Gaussian dressing matrices are `(d**2, d**2)` arrays. Vectorization column-stacks with `order="F"`. Prefer a density **matrix** for every high-level `rho0` input; a few low-level response functions additionally accept a vector as documented below.
- High-level calls first use an explicitly supplied non-`None` `rho0`, then `system.rho0`, and otherwise raise `ValueError`. There is no automatic stationary-state fallback. `rho_ss = system.stationary_state()` only returns a state: pass `rho0=rho_ss`, or supply that state when constructing the system.
- `L` generates `d rho/dt = L rho` in vectorized form. The Hamiltonian contribution is `-(1j/hbar)[H,rho]`, and the interaction is `+(1j/hbar)[mu,rho]`. Collapse operators include square-root rate factors, for example `sqrt(gamma)*C`.
- Resolvents use `G_eta(omega) = [(eta - 1j*omega)*I - L]^{-1}`, the causal transform of `exp(L*t)` weighted by `exp((1j*omega - eta)*t)`. For a stable physical generator, positive `eta` ensures convergence and regularizes undamped modes. It also adds broadening beyond physical dissipation. When it is only a convergence regulator, use a small positive value giving stable, adequately resolved results; reduce it while refining the grid and compare unnormalized complex responses. It need not be the smallest floating-point number. Frequency-domain routines do not generally validate positivity of `eta` or `T`, so supply physical values yourself.
- Third-order calls with `pathway`/`signature` use `"NR" -> (+1,-1,+1)` and `"R" -> (-1,+1,+1)`. An explicit three-sign `signature` overrides `pathway`. The signs select field components; only an RWA function also restricts the molecular interaction to `V_plus`/`V_minus`.
- Axes are signed and are never flipped automatically. For standard NR, use positive `omega1` and positive `omega3`; for standard R, use negative `omega1` and positive `omega3`. To display the usual positive R excitation coordinate, pass `-omega1_R` to plotting with the corresponding R array. Do not calculate R on a positive first-coherence axis merely by changing `pathway`.
- Every 2D spectrum has `S[j,i] = S(omega3[j], T, omega1[i])`. `finite_pulse_points` instead returns only pairs. The grid code caches separate frequency-dependent pieces; `response_grid` is the general scalar-loop alternative.
- `GaussianPulse` supplies carrier, width, peak amplitude and phase, not pulse arrival times. Frequency-domain functions use the ordered three-interaction construction and a separate waiting time `T`; arbitrary overlapping-field simulations are outside this API.
- The standard finite-pulse formulation assumes stationary preparation under the field-free Liouvillian. Accepting an explicitly supplied nonstationary `rho0` does not establish the formula's validity for an arbitrary prepared-state experiment. Check the scientific assumptions as well as the accepted array shapes.
- There are five third-order frequency-domain models but only four linear models. PULSUS 0.2.0 does **not** expose a `finite_pulse_linear_rwa_*` family. Do not infer an API name by symmetry.

## Parameter key

The exact signatures below show all required inputs (no default) and all optional defaults. Common names mean:

| Parameters | Meaning |
|---|---|
| `system` | `SpectroscopySystem` holding `L`, interactions, dipole, and optional state |
| `pulse`, `pulse1`, `pulse2`, `pulse3` | `GaussianPulse` objects; one object may be reused |
| `L_super` | dense field-free Liouvillian |
| `V`, `V_plus`, `V_minus` | full or signed interaction superoperators |
| `observable` | detection operator; high-level wrappers use `system.dipole` |
| `rho0` | initial density matrix, explicitly required in low-level responses |
| `omega` | scalar frequency for `*_signal`/`*_response`, 1D axis for linear spectra |
| `omega1`, `omega3` | scalar coordinates, 1D Cartesian axes, or paired arrays, according to the function |
| `T` | waiting-time interval between interactions 2 and 3 |
| `eta` | positive resolvent convergence/broadening rate |
| `omega_L`, `sigma`, `E0`, `phase` | one pulse's carrier, temporal width, peak amplitude and phase |
| `omega_L1`...`omega_L3`, `sigma1`...`sigma3`, `E01`...`E03`, `phase1`...`phase3` | low-level scalar pulse inputs; suffix indicates interaction number |
| `pathway`, `signature`, `sign` | standard alias, explicit three-sign sector, or one linear field sign |
| `times` | non-negative laboratory or delay time axis, as stated for that function |
| `t1`, `t3` | non-negative delay axes for a two-dimensional impulse-response grid |
| `pulse_time`, `pulse_times` | one laboratory arrival time, or ordered `(tau1,tau2,tau3)` |
| `F1`, `F2`, `F3` | explicitly supplied pulse-dressing matrices |

The implementation is dense NumPy/SciPy. Its matrix dimension grows as `d**2`; pulse dressing uses matrix exponentials. Begin with small Hilbert spaces and modest frequency grids.

## 1. Define a system and pulses

### `SpectroscopySystem`

```python
SpectroscopySystem(
    H,
    dipole,
    collapse_ops=None,
    rho0=None,
    hbar=1.0,
    dipole_plus=None,
    dipole_minus=None,
)
```

Build a dense Liouville-space system from a Hamiltonian, detection dipole, and optional Lindblad collapse operators.

**Required:** `H`, `dipole`.

**Returns:** A system object with `H`, `dipole`, `collapse_ops`, `rho0`, `hbar`, `d`, `L`, `V`, `dipole_plus`, `dipole_minus`, `V_plus`, and `V_minus` attributes.

`H` and `dipole` are `(d, d)` arrays. Collapse operators already include square-root rate factors. Supplying `dipole_plus` and `dipole_minus` independently constructs the two RWA interactions; the class does not derive or check that decomposition. The method `system.stationary_state()` returns a `(d, d)` state without changing `system.rho0`. Existing superoperators are built in the constructor: construct a new system after changing the physical operators.

[Source](../src/pulsus/system.py#L10).

### `GaussianPulse`

```python
GaussianPulse(omega_L: float, sigma: float, E0: float = 1.0, phase: float = 0.0)
```

Store one Gaussian pulse as an immutable dataclass.

**Required:** `omega_L`, `sigma`.

**Returns:** A pulse with `omega_L`, `sigma`, `E0`, and `phase` fields.

Use positive `sigma` in time units, carrier angular frequency `omega_L` in inverse time units, and phase in radians. `sigma` is the Gaussian field-envelope standard deviation; field-envelope FWHM is `2*sqrt(2*log(2))*sigma`, and intensity FWHM is `2*sqrt(log(2))*sigma`. No automatic width or unit validation is performed.

[Source](../src/pulsus/pulses.py#L6).

## 2. Calculate third-order frequency grids and selected points

### `finite_pulse_spectrum`

```python
finite_pulse_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with the complete interaction `V` on the Cartesian product of two frequency axes.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`; `S[j, i]` is the response at `(omega1[i], omega3[j])`.

Selecting field signs does not impose a molecular RWA: every interaction still uses the full `V`.

[Source](../src/pulsus/spectra.py#L186).

### `finite_pulse_rwa_spectrum`

```python
finite_pulse_rwa_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with signed RWA interactions `V_plus`/`V_minus` on the Cartesian product of two frequency axes.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`; `S[j, i]` is the response at `(omega1[i], omega3[j])`.

RWA calls require both `dipole_plus` and `dipole_minus` in the system.

[Source](../src/pulsus/spectra.py#L572).

### `short_pulse_rwa_spectrum`

```python
short_pulse_rwa_spectrum(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate scalar Gaussian spectra multiplying the impulsive RWA molecular response on the Cartesian product of two frequency axes.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`; `S[j, i]` is the response at `(omega1[i], omega3[j])`.

RWA calls require both `dipole_plus` and `dipole_minus` in the system. Pulse 2 uses `gaussian_spectrum(-omega1, ..., sign=s2)`; pulses 1 and 3 use `omega1` and `omega3`. This model retains finite spectral amplitudes but neglects molecular evolution during each pulse.

[Source](../src/pulsus/spectra.py#L766).

### `impulsive_spectrum`

```python
impulsive_spectrum(system, omega1, omega3, T, eta, rho0=None)
```

Evaluate the impulsive molecular response with the complete interaction `V` on the Cartesian product of two frequency axes.

**Required:** `system`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`; `S[j, i]` is the response at `(omega1[i], omega3[j])`.

No `pathway`, `signature`, or pulse arguments are accepted. The full molecular response is evaluated at the frequencies supplied.

[Source](../src/pulsus/spectra.py#L23).

### `impulsive_rwa_spectrum`

```python
impulsive_rwa_spectrum(
    system,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate the impulsive molecular response with signed RWA interactions on the Cartesian product of two frequency axes.

**Required:** `system`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`; `S[j, i]` is the response at `(omega1[i], omega3[j])`.

RWA calls require both `dipole_plus` and `dipole_minus` in the system.

[Source](../src/pulsus/spectra.py#L866).

### `finite_pulse_points`

```python
finite_pulse_points(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate the full-interaction finite-pulse model only at paired coordinates `(omega1[k], omega3[k])`, caching repeated frequencies.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex one-dimensional array of shape `(N,)`, in the original pair order.

Both frequency inputs must be one-dimensional and have identical shape. These are pairs, not a Cartesian product. There is no separate RWA paired-point function in 0.2.0.

[Source](../src/pulsus/spectra.py#L365).

### `response_grid`

```python
response_grid(response_function, omega1, omega3, **kwargs)
```

Evaluate an arbitrary scalar response callable at every pair on two one-dimensional frequency axes.

**Required:** `response_function`, `omega1`, `omega3`.

**Returns:** A complex array of shape `(len(omega3), len(omega1))`.

`response_function` must accept `omega1` and `omega3` keyword arguments and return a complex scalar; other required inputs go in `**kwargs`. This generic nested-loop evaluator does not use the caching in the specialized spectrum functions.

[Source](../src/pulsus/spectra.py#L138).

## 3. Calculate a third-order value at one coordinate

### `finite_pulse_signal`

```python
finite_pulse_signal(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with the complete interaction `V` at one frequency pair.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

Selecting field signs does not impose a molecular RWA: every interaction still uses the full `V`.

[Source](../src/pulsus/response.py#L629).

### `finite_pulse_rwa_signal`

```python
finite_pulse_rwa_signal(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with signed RWA interactions `V_plus`/`V_minus` at one frequency pair.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

RWA calls require both `dipole_plus` and `dipole_minus` in the system.

[Source](../src/pulsus/response.py#L681).

### `short_pulse_rwa_signal`

```python
short_pulse_rwa_signal(
    system,
    pulse1,
    pulse2,
    pulse3,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate scalar Gaussian spectra multiplying the impulsive RWA molecular response at one frequency pair.

**Required:** `system`, `pulse1`, `pulse2`, `pulse3`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

RWA calls require both `dipole_plus` and `dipole_minus` in the system. Pulse 2 uses `gaussian_spectrum(-omega1, ..., sign=s2)`; pulses 1 and 3 use `omega1` and `omega3`. This model retains finite spectral amplitudes but neglects molecular evolution during each pulse.

[Source](../src/pulsus/response.py#L784).

### `impulsive_signal`

```python
impulsive_signal(system, omega1, omega3, T, eta, rho0=None)
```

Evaluate the impulsive molecular response with the complete interaction `V` at one frequency pair.

**Required:** `system`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

No `pathway`, `signature`, or pulse arguments are accepted. The full molecular response is evaluated at the frequencies supplied.

[Source](../src/pulsus/response.py#L842).

### `impulsive_rwa_signal`

```python
impulsive_rwa_signal(
    system,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    rho0=None,
    signature=None,
)
```

Evaluate the impulsive molecular response with signed RWA interactions at one frequency pair.

**Required:** `system`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

RWA calls require both `dipole_plus` and `dipole_minus` in the system.

[Source](../src/pulsus/response.py#L740).

## 4. Calculate and plot linear spectra

### `finite_pulse_linear_spectrum`

```python
finite_pulse_linear_spectrum(system, pulse, omega, eta, sign=1, rho0=None)
```

Evaluate the full-interaction finite-Gaussian response using a system object.

**Required:** `system`, `pulse`, `omega`, `eta`.

**Returns:** A complex array of shape `(len(omega),)`.

`omega` must be one-dimensional. The `sign` selects a field-frequency component while the molecular interaction remains the complete `V`; this is not molecular RWA.

[Source](../src/pulsus/linear.py#L411).

### `impulsive_linear_spectrum`

```python
impulsive_linear_spectrum(system, omega, eta, rho0=None)
```

Evaluate the full-interaction impulsive molecular response using a system object.

**Required:** `system`, `omega`, `eta`.

**Returns:** A complex array of shape `(len(omega),)`.

`omega` must be one-dimensional. 

[Source](../src/pulsus/linear.py#L378).

### `impulsive_linear_rwa_spectrum`

```python
impulsive_linear_rwa_spectrum(system, omega, eta, sign=1, rho0=None)
```

Evaluate the impulsive RWA molecular response using a system object.

**Required:** `system`, `omega`, `eta`.

**Returns:** A complex array of shape `(len(omega),)`.

`omega` must be one-dimensional. RWA wrappers require both system dipole components. `sign=+1` selects `V_plus`; `sign=-1` selects `V_minus`.

[Source](../src/pulsus/linear.py#L449).

### `short_pulse_linear_rwa_spectrum`

```python
short_pulse_linear_rwa_spectrum(system, pulse, omega, eta, sign=1, rho0=None)
```

Evaluate the scalar-Gaussian-dressed impulsive RWA response using a system object.

**Required:** `system`, `pulse`, `omega`, `eta`.

**Returns:** A complex array of shape `(len(omega),)`.

`omega` must be one-dimensional. RWA wrappers require both system dipole components. The result is `gaussian_spectrum(omega, ...)` times the impulsive RWA molecular response.

[Source](../src/pulsus/linear.py#L484).

### `finite_pulse_linear_signal`

```python
finite_pulse_linear_signal(system, pulse, omega, eta, sign=1, rho0=None)
```

Evaluate the full-interaction finite-Gaussian response using a system object.

**Required:** `system`, `pulse`, `omega`, `eta`.

**Returns:** A complex scalar.

The `sign` selects a field-frequency component while the molecular interaction remains the complete `V`; this is not molecular RWA.

[Source](../src/pulsus/linear.py#L251).

### `impulsive_linear_signal`

```python
impulsive_linear_signal(system, omega, eta, rho0=None)
```

Evaluate the full-interaction impulsive molecular response using a system object.

**Required:** `system`, `omega`, `eta`.

**Returns:** A complex scalar.



[Source](../src/pulsus/linear.py#L223).

### `impulsive_linear_rwa_signal`

```python
impulsive_linear_rwa_signal(system, omega, eta, sign=1, rho0=None)
```

Evaluate the impulsive RWA molecular response using a system object.

**Required:** `system`, `omega`, `eta`.

**Returns:** A complex scalar.

RWA wrappers require both system dipole components. `sign=+1` selects `V_plus`; `sign=-1` selects `V_minus`.

[Source](../src/pulsus/linear.py#L289).

### `short_pulse_linear_rwa_signal`

```python
short_pulse_linear_rwa_signal(system, pulse, omega, eta, sign=1, rho0=None)
```

Evaluate the scalar-Gaussian-dressed impulsive RWA response using a system object.

**Required:** `system`, `pulse`, `omega`, `eta`.

**Returns:** A complex scalar.

RWA wrappers require both system dipole components. The result is `gaussian_spectrum(omega, ...)` times the impulsive RWA molecular response.

[Source](../src/pulsus/linear.py#L329).

## 5. Display a third-order spectrum

### `plot_spectrum`

```python
plot_spectrum(
    omega1,
    omega3,
    spectrum,
    component='real',
    normalize=False,
    levels=31,
    ax=None,
    colorbar=True,
    colorbar_label=None,
    title=None,
    xlabel='$\\omega_1$',
    ylabel='$\\omega_3$',
)
```

Display one component of a complex 2D grid using Matplotlib filled contours.

**Required:** `omega1`, `omega3`, `spectrum`.

**Returns:** `(fig, ax, contour)`: a Figure, Axes, and QuadContourSet.

`spectrum` must have shape `(len(omega3), len(omega1))`. `component` is `"real"`, `"imag"`, or `"abs"`. `normalize=True` divides only the displayed component by its own maximum absolute value (zero remains zero); it does not alter the calculated complex spectrum. Signed components use a symmetric diverging scale. This function neither phases nor constructs an absorptive spectrum, and does not save or show the figure automatically.

[Source](../src/pulsus/plotting.py#L5).

## 6. Inspect impulsive time-domain responses and traces

### `impulsive_linear_time_signal`

```python
impulsive_linear_time_signal(system, times, rho0=None)
```

Evaluate the full-interaction impulse response `<mu| exp(L*t) V |rho0>` at delay times after the interaction.

**Required:** `system`, `times`.

**Returns:** A complex array of shape `(len(times),)`.

`times` must be a one-dimensional non-negative axis. No pulse envelope, field amplitude, or broadening `eta` is applied.

[Source](../src/pulsus/time_domain.py#L48).

### `impulsive_linear_rwa_time_signal`

```python
impulsive_linear_rwa_time_signal(system, times, sign=1, rho0=None)
```

Evaluate one RWA impulse-response sector at non-negative delay times.

**Required:** `system`, `times`.

**Returns:** A complex array of shape `(len(times),)`.

Requires both RWA dipole components; `sign` must be +1 or -1. No Gaussian envelope or `eta` is applied.

[Source](../src/pulsus/time_domain.py#L92).

### `impulsive_third_order_time_signal`

```python
impulsive_third_order_time_signal(system, t1, t3, T, rho0=None)
```

Evaluate the full-interaction impulse-response grid `<mu| exp(L*t3) V exp(L*T) V exp(L*t1) V |rho0>`.

**Required:** `system`, `t1`, `t3`, `T`.

**Returns:** A complex array of shape `(len(t3), len(t1))`.

`t1` and `t3` must be one-dimensional non-negative axes; `T` must be non-negative. No sign sector, Gaussian envelope, or `eta` is applied.

[Source](../src/pulsus/time_domain.py#L226).

### `impulsive_third_order_rwa_time_signal`

```python
impulsive_third_order_rwa_time_signal(
    system,
    t1,
    t3,
    T,
    pathway='NR',
    signature=None,
    rho0=None,
)
```

Evaluate one signed RWA third-order impulse-response grid.

**Required:** `system`, `t1`, `t3`, `T`.

**Returns:** A complex array of shape `(len(t3), len(t1))`.

Same time-axis conventions as the full-interaction grid. Requires both RWA dipole components. Note the positional order of optional arguments: `pathway`, `signature`, then `rho0`. Prefer keywords.

[Source](../src/pulsus/time_domain.py#L255).

### `impulsive_linear_polarization`

```python
impulsive_linear_polarization(system, times, pulse_time=0.0, sign=None, rho0=None)
```

Evaluate a laboratory-time first-order impulsive trace, zero before `pulse_time`.

**Required:** `system`, `times`.

**Returns:** A complex array of shape `(len(times),)`.

`sign=None` uses full `V`; +1 or -1 selects RWA and requires both dipole components. Times and pulse time must be non-negative. The supplied state is used at the interaction; no field-free evolution from laboratory time zero to `pulse_time` is performed.

[Source](../src/pulsus/time_domain.py#L304).

### `impulsive_third_order_polarization`

```python
impulsive_third_order_polarization(
    system,
    times,
    pulse_times,
    pathway=None,
    signature=None,
    rho0=None,
)
```

Evaluate a laboratory-time third-order impulsive trace for one ordered three-pulse sequence, zero before the third pulse.

**Required:** `system`, `times`, `pulse_times`.

**Returns:** A complex array of shape `(len(times),)`.

`pulse_times=(tau1, tau2, tau3)` must be non-negative with `tau1 <= tau2 <= tau3`. `pathway=None` and `signature=None` use full `V`; supplying either selects RWA. The state is interpreted at the first interaction; for a nonstationary state prepared earlier, evolve it yourself to `tau1`. The trace does not sum lower-order stages.

[Source](../src/pulsus/time_domain.py#L383).

### `impulsive_pathway_polarization`

```python
impulsive_pathway_polarization(
    system,
    times,
    pulse_times,
    pathway='NR',
    signature=None,
    rho0=None,
    return_components=False,
)
```

Evaluate the piecewise laboratory-time history of one ordered impulsive RWA sign sector.

**Required:** `system`, `times`, `pulse_times`.

**Returns:** A complex array of shape `(len(times),)`; with `return_components=True`, returns `(polarization, components)` where the dict keys are `"first_order"`, `"second_order"`, and `"third_order"`, each a complex array of that shape.

Requires both RWA dipole components. Pulse ordering and state timing match the third-order trace. Stages occupy `[tau1,tau2)`, `[tau2,tau3)`, and `[tau3,infinity)`, respectively, and are zero outside their intervals. The piecewise sum describes different perturbative orders in different intervals; it is not a complete field-driven polarization simulation. In a simple electronic model the detected second-order stage can be zero.

[Source](../src/pulsus/time_domain.py#L546).

## 7. Select field signs

### `pathway_signs`

```python
pathway_signs(pathway)
```

Translate a standard R/NR alias to its three field-frequency signs.

**Required:** `pathway`.

**Returns:** `(+1, -1, +1)` for `"NR"`; `(-1, +1, +1)` for `"R"`.

Aliases are case insensitive. Other strings raise `ValueError`. This helper does not choose or transform frequency axes.

[Source](../src/pulsus/response.py#L282).

### `resolve_signature`

```python
resolve_signature(pathway='NR', signature=None)
```

Resolve an explicit field-sign tuple or fall back to a standard pathway alias.

**Required:** none.

**Returns:** A three-element tuple `(s1, s2, s3)`.

An explicit `signature` overrides `pathway`; each entry must be +1 or -1 and length must be three. There are eight possible signatures. An invalid `pathway` is not inspected if a valid signature is supplied. A sign sector contains the molecular contributions selected by the interactions; it is not a single microscopic double-sided Feynman diagram.

[Source](../src/pulsus/response.py#L307).

## 8. Build Liouville-space operators and resolvents

### `vec`

```python
vec(rho)
```

Column-stack a matrix, using NumPy reshape order `"F"`.

**Required:** `rho`.

**Returns:** A complex one-dimensional vector; for a `(d, d)` density matrix its length is `d**2`.

The helper converts its input to a complex NumPy array; it does not check that a matrix is a physical density operator.

[Source](../src/pulsus/liouville.py#L4).

### `unvec`

```python
unvec(rho_vec, d)
```

Restore a column-stacked vector to a square matrix.

**Required:** `rho_vec`, `d`.

**Returns:** A complex `(d, d)` array.

The vector must contain exactly `d**2` entries.

[Source](../src/pulsus/liouville.py#L22).

### `hamiltonian_liouvillian`

```python
hamiltonian_liouvillian(H, hbar=1.0)
```

Build the Hamiltonian generator `rho -> -(1j/hbar)*(H@rho - rho@H)`.

**Required:** `H`.

**Returns:** A dense complex `(d**2, d**2)` array.

Uses column stacking; `H` is `(d, d)` and must use energy units consistent with `hbar` and the chosen time unit.

[Source](../src/pulsus/liouville.py#L42).

### `interaction_superoperator`

```python
interaction_superoperator(mu, hbar=1.0)
```

Build the interaction `rho -> (1j/hbar)*(mu@rho - rho@mu)`.

**Required:** `mu`.

**Returns:** A dense complex `(d**2, d**2)` array.

The same constructor can be used for the full dipole or either RWA component; `mu` is a `(d, d)` operator.

[Source](../src/pulsus/liouville.py#L78).

### `dissipator_super`

```python
dissipator_super(L, I=None)
```

Build one Lindblad dissipator `C@rho@C.conj().T - 0.5*{C.conj().T@C,rho}`.

**Required:** `L`.

**Returns:** A dense complex `(d**2, d**2)` array.

The parameter named `L` here is the `(d, d)` collapse/jump operator, not the Liouvillian `L_super`. Include the square root of the rate in it. Optional `I` is the Hilbert-space identity `(d, d)`.

[Source](../src/pulsus/liouville.py#L114).

### `liouvillian`

```python
liouvillian(H, collapse_ops=None, hbar=1.0)
```

Combine Hamiltonian evolution with the dissipators of all collapse operators.

**Required:** `H`.

**Returns:** A dense complex `(d**2, d**2)` array.

`collapse_ops=None` gives Hamiltonian-only evolution. Rates must already be included in the collapse operators; no extra rate argument is accepted.

[Source](../src/pulsus/liouville.py#L151).

### `stationary_state`

```python
stationary_state(L_super, d)
```

Select the eigenvector of `L_super` whose eigenvalue is closest to zero, restore its matrix shape, normalize its trace, and symmetrize it.

**Required:** `L_super`, `d`.

**Returns:** A complex `(d, d)` density-matrix candidate with unit trace.

Pass the Hilbert-space dimension `d`, not `d**2`. The routine does not guarantee a unique or positive stationary state and does not verify the residual. Check `norm(L_super @ vec(rho_ss))`, Hermiticity, trace, and eigenvalues when needed, especially for degenerate stationary subspaces.

[Source](../src/pulsus/liouville.py#L181).

### `observable_bra`

```python
observable_bra(observable)
```

Build the row coefficients that contract a vectorized state to `Tr(observable @ rho)`.

**Required:** `observable`.

**Returns:** A complex one-dimensional array of length `d**2`.

The result is `vec(observable.conj().T).conj()`. Contract it with `@`, using the package column-stacking convention.

[Source](../src/pulsus/response.py#L156).

### `resolvent_action`

```python
resolvent_action(L_super, omega, rhs, eta)
```

Solve `((eta - 1j*omega)*I - L_super) @ x = rhs` without constructing an inverse.

**Required:** `L_super`, `omega`, `rhs`, `eta`.

**Returns:** A complex array with the same vector or matrix right-hand-side shape as `rhs`.

`L_super` is `(d**2, d**2)`; `rhs` can be length `d**2` or `(d**2, m)`. Choose positive `eta` for convergence/broadening. A singular coefficient matrix can raise a SciPy linear-algebra error.

[Source](../src/pulsus/resolvents.py#L5).

### `resolvent_matrix`

```python
resolvent_matrix(L_super, omega, eta)
```

Construct the complete Liouvillian resolvent by solving against an identity matrix.

**Required:** `L_super`, `omega`, `eta`.

**Returns:** A dense complex `(d**2, d**2)` matrix.

For applying a resolvent to a state, prefer `resolvent_action`. Same signed-frequency and positive-eta conventions apply.

[Source](../src/pulsus/resolvents.py#L46).

## 9. Evaluate Gaussian field and dressing functions

### `gaussian_prefactor`

```python
gaussian_prefactor(sigma, E0=1.0, phase=0.0, sign=1)
```

Evaluate the prefactor of one signed Gaussian field-frequency component.

**Required:** `sigma`.

**Returns:** A scalar `0.5*E0*sqrt(2*pi)*sigma*exp(1j*sign*phase)`.

`sign` must be +1 or -1. At fixed peak amplitude `E0`, this prefactor scales with `sigma`; varying the width does not keep the pulse area fixed.

[Source](../src/pulsus/pulses.py#L27).

### `gaussian_spectrum`

```python
gaussian_spectrum(omega, omega_L, sigma, E0=1.0, phase=0.0, sign=1)
```

Evaluate one signed Gaussian field spectrum centered on `sign*omega_L`.

**Required:** `omega`, `omega_L`, `sigma`.

**Returns:** A scalar, or a NumPy array matching a NumPy-array `omega` input.

The formula is `gaussian_prefactor(...) * exp(-0.5*sigma**2*(omega - sign*omega_L)**2)`. It is one field-frequency component, not the sum of both signs. Supply a NumPy array for vector input.

[Source](../src/pulsus/pulses.py#L61).

### `pulse_dressing_1`

```python
pulse_dressing_1(L_super, omega1, omega_L, sigma, E0=1.0, phase=0.0, sign=1)
```

Construct the Gaussian Liouvillian dressing matrix for pulse 1.

**Required:** `L_super`, `omega1`, `omega_L`, `sigma`.

**Returns:** A dense complex `(d**2, d**2)` matrix.

Uses `F = gaussian_prefactor(...) * expm(0.5*sigma**2*(L_super + 1j*delta*I) @ (L_super + 1j*delta*I))`, with `delta = omega1 - sign*omega_L`. The squared matrix inside the exponential is a matrix product. 

[Source](../src/pulsus/pulses.py#L127).

### `pulse_dressing_2`

```python
pulse_dressing_2(L_super, omega1, omega_L, sigma, E0=1.0, phase=0.0, sign=1)
```

Construct the Gaussian Liouvillian dressing matrix for pulse 2.

**Required:** `L_super`, `omega1`, `omega_L`, `sigma`.

**Returns:** A dense complex `(d**2, d**2)` matrix.

Uses `F = gaussian_prefactor(...) * expm(0.5*sigma**2*(L_super + 1j*delta*I) @ (L_super + 1j*delta*I))`, with `delta = omega1 + sign*omega_L`. The squared matrix inside the exponential is a matrix product. Pulse 2 depends on `omega1`, not `omega3`.

[Source](../src/pulsus/pulses.py#L156).

### `pulse_dressing_3`

```python
pulse_dressing_3(L_super, omega3, omega_L, sigma, E0=1.0, phase=0.0, sign=1)
```

Construct the Gaussian Liouvillian dressing matrix for pulse 3.

**Required:** `L_super`, `omega3`, `omega_L`, `sigma`.

**Returns:** A dense complex `(d**2, d**2)` matrix.

Uses `F = gaussian_prefactor(...) * expm(0.5*sigma**2*(L_super + 1j*delta*I) @ (L_super + 1j*delta*I))`, with `delta = omega3 - sign*omega_L`. The squared matrix inside the exponential is a matrix product. 

[Source](../src/pulsus/pulses.py#L185).

## 10. Evaluate low-level third-order responses

### `finite_pulse_response`

```python
finite_pulse_response(
    L_super,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma1,
    sigma2,
    sigma3,
    pathway='NR',
    E01=1.0,
    E02=1.0,
    E03=1.0,
    phase1=0.0,
    phase2=0.0,
    phase3=0.0,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with the complete interaction `V` at one frequency pair.

**Required:** `L_super`, `V`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`, `omega_L1`, `omega_L2`, `omega_L3`, `sigma1`, `sigma2`, `sigma3`.

**Returns:** A complex scalar.

Inputs are explicit Liouville-space operators rather than a system/pulse object. Pass `rho0` as a `(d, d)` density matrix. The implementation vectorizes it directly. Pulse parameters are separate scalar inputs, with optional amplitudes and phases.

[Source](../src/pulsus/response.py#L356).

### `finite_pulse_rwa_response`

```python
finite_pulse_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma1,
    sigma2,
    sigma3,
    pathway='NR',
    E01=1.0,
    E02=1.0,
    E03=1.0,
    phase1=0.0,
    phase2=0.0,
    phase3=0.0,
    signature=None,
)
```

Evaluate finite Gaussian operator dressing with signed RWA interactions `V_plus`/`V_minus` at one frequency pair.

**Required:** `L_super`, `V_plus`, `V_minus`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`, `omega_L1`, `omega_L2`, `omega_L3`, `sigma1`, `sigma2`, `sigma3`.

**Returns:** A complex scalar.

Inputs are explicit Liouville-space operators rather than a system/pulse object. Pass `rho0` as a `(d, d)` density matrix. The implementation vectorizes it directly. Pulse parameters are separate scalar inputs, with optional amplitudes and phases.

[Source](../src/pulsus/response.py#L455).

### `short_pulse_rwa_response`

```python
short_pulse_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    omega_L1,
    omega_L2,
    omega_L3,
    sigma1,
    sigma2,
    sigma3,
    pathway='NR',
    E01=1.0,
    E02=1.0,
    E03=1.0,
    phase1=0.0,
    phase2=0.0,
    phase3=0.0,
    signature=None,
)
```

Evaluate scalar Gaussian spectra multiplying the impulsive RWA molecular response at one frequency pair.

**Required:** `L_super`, `V_plus`, `V_minus`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`, `omega_L1`, `omega_L2`, `omega_L3`, `sigma1`, `sigma2`, `sigma3`.

**Returns:** A complex scalar.

Inputs are explicit Liouville-space operators rather than a system/pulse object. Accepts a density matrix or column-stacked vector via the impulsive RWA molecular response. Pulse 2 is evaluated at `-omega1`.

[Source](../src/pulsus/response.py#L14).

### `impulsive_response`

```python
impulsive_response(L_super, V, observable, rho0, omega1, omega3, T, eta)
```

Evaluate the impulsive molecular response with the complete interaction `V` at one frequency pair.

**Required:** `L_super`, `V`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

Inputs are explicit Liouville-space operators rather than a system/pulse object. Accepts a density matrix or a column-stacked vector. No sign sector is selected.

[Source](../src/pulsus/response.py#L563).

### `impulsive_rwa_response`

```python
impulsive_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    pathway='NR',
    signature=None,
)
```

Evaluate the impulsive molecular response with signed RWA interactions at one frequency pair.

**Required:** `L_super`, `V_plus`, `V_minus`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`.

**Returns:** A complex scalar.

Inputs are explicit Liouville-space operators rather than a system/pulse object. Accepts a density matrix or a column-stacked vector; the signs select `V_plus` or `V_minus` at each interaction.

[Source](../src/pulsus/response.py#L91).

### `third_order_response`

```python
third_order_response(
    L_super,
    V,
    observable,
    rho0,
    omega1,
    omega3,
    T,
    eta,
    F1,
    F2,
    F3,
)
```

Evaluate the full-interaction finite-pulse chain using three precomputed dressing matrices.

**Required:** `L_super`, `V`, `observable`, `rho0`, `omega1`, `omega3`, `T`, `eta`, `F1`, `F2`, `F3`.

**Returns:** A complex scalar.

`F1`, `F2`, and `F3` must be `(d**2, d**2)` matrices for the requested frequencies. Accepts a density matrix or column-stacked vector. The operator order is `<observable| G3 V exp(L*T) F3 F2 V G1 V F1 |rho0>`; `Gk = [(eta - 1j*omegak)I - L]**(-1)` here denotes an inverse, not elementwise power.

[Source](../src/pulsus/response.py#L179).

## 11. Evaluate low-level linear responses

### `finite_pulse_linear_response`

```python
finite_pulse_linear_response(
    L_super,
    V,
    observable,
    rho0,
    omega,
    eta,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
)
```

Evaluate the full-interaction finite-Gaussian response using explicit operators.

**Required:** `L_super`, `V`, `observable`, `rho0`, `omega`, `eta`, `omega_L`, `sigma`.

**Returns:** A complex scalar.

The `sign` selects a field-frequency component while the molecular interaction remains the complete `V`; this is not molecular RWA.

[Source](../src/pulsus/linear.py#L73).

### `impulsive_linear_response`

```python
impulsive_linear_response(L_super, V, observable, rho0, omega, eta)
```

Evaluate the full-interaction impulsive molecular response using explicit operators.

**Required:** `L_super`, `V`, `observable`, `rho0`, `omega`, `eta`.

**Returns:** A complex scalar.



[Source](../src/pulsus/linear.py#L38).

### `impulsive_linear_rwa_response`

```python
impulsive_linear_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega,
    eta,
    sign=1,
)
```

Evaluate the impulsive RWA molecular response using explicit operators.

**Required:** `L_super`, `V_plus`, `V_minus`, `observable`, `rho0`, `omega`, `eta`.

**Returns:** A complex scalar.

`sign=+1` selects `V_plus`; `sign=-1` selects `V_minus`.

[Source](../src/pulsus/linear.py#L135).

### `short_pulse_linear_rwa_response`

```python
short_pulse_linear_rwa_response(
    L_super,
    V_plus,
    V_minus,
    observable,
    rho0,
    omega,
    eta,
    omega_L,
    sigma,
    E0=1.0,
    phase=0.0,
    sign=1,
)
```

Evaluate the scalar-Gaussian-dressed impulsive RWA response using explicit operators.

**Required:** `L_super`, `V_plus`, `V_minus`, `observable`, `rho0`, `omega`, `eta`, `omega_L`, `sigma`.

**Returns:** A complex scalar.

The result is `gaussian_spectrum(omega, ...)` times the impulsive RWA molecular response.

[Source](../src/pulsus/linear.py#L173).

## Completeness and scope

The 57 entries above include both classes, all high-level frequency and time routines, plotting, every root-imported low-level response, and all root-imported Liouville, resolvent, Gaussian, and signature helpers. Private underscore-prefixed implementation helpers are not part of this reference. Code in `examples/` and notebooks may define additional local helpers; those are not package-root API.

A source update should trigger a documentation review of `src/pulsus/__init__.py`, the function signatures, shapes, and examples. The API count here is a source audit, not a claim that all scientific regimes have been independently validated by this manual.
