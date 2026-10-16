# Advanced numerical notes

These optional notes collect material moved out of the introductory user manual. For the recorded repository checks, see [VALIDATION.md](VALIDATION.md).

## Comparing impulsive time and frequency results

The resolvent is the one-sided transform

```text
G_eta(omega) = integral_0^infinity exp[(i*omega - eta)*t] exp[L*t] dt
             = [(eta - i*omega)*I - L]^-1.
```

For a stable physical Lindblad generator and positive `eta`, this full operator integral converges. The generator has stationary modes and may contain undamped coherences; positive `eta` also regularizes their resolvent poles. Physical dissipation can already make a particular response channel decay, but convergence of that channel is distinct from invertibility of the full matrix solved by the implementation.

For one isolated mode with eigenvalue `-Gamma - i*Omega`, the resolvent denominator is `Gamma + eta - i*(omega - Omega)`. Thus `eta` adds to its Lorentzian half-width. Multipeak interference and different plotted components need not follow one simple apparent-linewidth rule.

When `eta` is only a convergence regulator, reduce it while resolving the resulting peaks and check that unnormalized complex spectra are stable. The aim is small added broadening with stable linear solves, not the smallest machine-representable number. With no physical linewidth, the `eta -> 0+` limit can be singular or distributional, so finite peak-height convergence is not guaranteed; specify physical damping or report the broadening/resolution chosen.

For a third-order comparison, transform both `t1` and `t3` with `exp[i*(omega1*t1 + omega3*t3)] * exp[-eta*(t1+t3)]` at fixed waiting time `T`. Do not add `exp[-eta*T]` when reproducing the frequency API. Finite integration windows, timestep, quadrature and FFT sign conventions introduce additional numerical errors. The direct frequency functions avoid this time-sampling step.

## Checking a new model

- Verify compatible operator shapes, Hermitian Hamiltonian and dipole, and a normalized, Hermitian, positive semidefinite initial density matrix.
- For stationary preparation, inspect `norm(system.L @ pulsus.vec(rho0))` relative to the scale of `system.L`. No single absolute tolerance applies to every unit system.
- Check finite spectral values and use unequal excitation and detection grid lengths to reveal accidental axis swaps.
- Compare a few grid entries with the corresponding scalar or paired-coordinate calculation, using identical inputs.
- Refine spectral spacing around narrow features and vary `eta` independently of the physical collapse operators.

The stationary-state solver selects the eigenvector whose eigenvalue is closest to zero, normalizes its trace and symmetrizes it. It does not prove uniqueness or enforce positive semidefiniteness. A purely Hamiltonian model or a model with disconnected sectors can have multiple stationary states.

## Comparing approximations

For complex arrays on identical coordinates, one useful error is `norm(S_test - S_ref) / norm(S_ref)`, using a nonzero reference norm. Report the spectral region and whether a global complex scale was fitted. Normalizing each response separately can conceal amplitude differences; fitting a complex scale can conceal both amplitude and phase differences.

For finite pulses, reassess the separated-pulse assumption when changing widths or delays. For truncated vibronic systems, increase the vibrational basis cutoff as well as refining the frequency grid. Keep the same initial-state prescription and dipole definition across model comparisons. Holding peak amplitude fixed and holding envelope area fixed are different pulse comparisons.
