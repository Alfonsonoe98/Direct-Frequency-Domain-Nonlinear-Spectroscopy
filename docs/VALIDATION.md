# Documentation validation record

Manual version 0.3, targeting PULSUS 0.2.0. Prepared October 15, 2026.

## Scope of the checks

- The existing repository suite passed: **67 tests**.
- All five scripts in `examples/manual/` completed. They display figures with `plt.show()`; the separate maintainer renderer regenerated all five PDF figures.
- Executable Python blocks in the usage chapters were checked in their stated tutorial context. Signature-only blocks are reference material; physical-unit templates were supplied with consistent example inputs for the check.
- The complete Python block in `docs/quickstart.md` ran by importing the bundled `pulsus_spec-0.2.0-py3-none-any.whl`, independently of the repository package helpers. The wheel's package modules match the supplied source.
- The selected-point example agreed with the corresponding Cartesian grid entries to a maximum absolute difference of about `2.3e-13`.
- The API Markdown reference contains all 57 root-imported names and source-matched signatures. Two of these are omitted from the release's `__all__`; the documentation records this without changing library code.
- LaTeX compilation and rendered-page review were used to check the PDF.

The package-only check used the local wheel with the existing environment's dependencies. It was not a fresh download/install from PyPI, and no claim of a fully isolated dependency installation is made. Passing these checks does not establish physical convergence for a different system or approximation regime.

## Environment used

Python 3.11.4

| Package | Version |
|---|---|
| pulsus-spec | 0.2.0 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| matplotlib | 3.11.1 |
| qutip | 5.3.1 |
| pytest | 9.1.1 |

## Repeat the checks

From the repository root in an active editable development environment:

```bash
python -m pytest -q
python examples/manual/01_first_spectrum.py
python examples/manual/02_linear_models.py
python examples/manual/03_third_order_models.py
python examples/manual/04_selected_points.py
python examples/manual/05_time_domain.py
python docs/user_manual/render_figures.py
cd docs/user_manual
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```
