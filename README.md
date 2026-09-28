# Direct Frequency-Domain Nonlinear Spectroscopy

Development and numerical validation of a direct frequency-domain
formulation of nonlinear spectroscopy for finite-dimensional open
quantum systems.

The project uses Liouvillian resolvents for direct evaluation of
frequency-domain response functions and finite-width Gaussian pulse
operators for treating molecular evolution during the laser pulses.

## Current status

- [x] Liouville-space formulation
- [x] Direct third-order frequency-domain response
- [x] Short and separated pulse formulation
- [x] Finite-width Gaussian pulse operators
- [x] Coupled two-level thermal model defined
- [ ] Coupled-dimer numerical implementation
- [ ] Linear-response validation
- [ ] Third-order time-domain benchmark
- [ ] Rephasing/nonrephasing benchmark
- [ ] Finite-width pulse benchmark
- [ ] Pulse-width sweep

## Repository structure

- `main.tex` — project technical note
- `references.bib` — bibliography database
- `literature_notes.md` — annotated literature notes
- `notebooks/` — numerical validation notebooks
- `src/` — reusable Python implementation
- `tests/` — numerical and consistency tests
- `figures/` — manuscript figures
- `results/` — compact benchmark results