# Documentation and package validation record

Manual version 0.5, targeting PULSUS 0.3.0. Prepared October 17, 2026.

## Checks completed

- The complete source-checkout suite passed: **108 tests**.
- A wheel and source distribution were built in isolation. The wheel and its dependencies were installed in a fresh virtual environment; the complete suite also passed against the installed wheel. The existing test-only `dimer_model.py` fixture was supplied separately, without placing the source package on the import path.
- The 46 plotting cases include component values, normalization and zero responses, input preservation, axes reuse, optional PNG/PDF export, actual raster DPI, continued editing and saving, and compatibility with existing 2D positional calls.
- All five scripts in `examples/manual/` ran against the installed wheel from an empty working folder. They displayed figures through the usual `plt.show()` call and wrote no files by default.
- All quick-start Python blocks, both complete README examples, and the API overlay/export example were executed in their documented context.
- The maintainer renderer regenerated the linear-comparison figure without point markers. The other four manual figures were preserved because their calculations and presentation did not change.
- The Markdown reference covers all 58 root-imported names and matching plotting signatures. The two pre-existing omissions from `__all__` remain documented; this release adds the new plotter without unrelated export changes.
- All 109 checked Markdown links resolve to existing local or repository targets. The shortened README preserves the original 24 display equations and scientific benchmarks in [method notes](method_notes.md). LaTeX compilation completed without undefined references or overfull boxes, and all 26 rendered pages were visually reviewed.

These checks verify implementation and documentation consistency. They do not establish physical convergence for arbitrary systems or approximation regimes. They were performed on the prepared 0.3.0 distributions before publication; verification of the published PyPI artifacts is a separate release step.

## Published release verification

On October 17, 2026, the GitHub publication workflow completed successfully for tag `v0.3.0`, commit `e1989ad3b262ef2e46ac968682bc4ecfbd0c358c`. [PyPI](https://pypi.org/project/pulsus-spec/0.3.0/) provides the wheel and source distribution; [Zenodo](https://doi.org/10.5281/zenodo.23432503) archives the release.

A new isolated environment installed `pulsus-spec==0.3.0` from PyPI, rather than the source checkout or locally built wheel. All ten installed Python source files match the tagged release. A two-level-system calculation produced finite linear and 2D responses. Both plotters saved nothing by default, returned editable figures, and exported valid PNG/PDF files through their optional arguments and through later `fig.savefig(...)` calls. The linear plotter had no markers by default.

The version DOI is `10.5281/zenodo.23432503`. The all-versions concept DOI is `10.5281/zenodo.22755139`; Zenodo's public metadata identifies the archived version as `v0.3.0`, published October 17, 2026.

## Environment

Python 3.11.4

| Package | Version |
|---|---|
| pulsus-spec | 0.3.0 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| matplotlib | 3.11.2 |
| qutip | 5.3.1 |
| pytest | 9.1.1 |

## Repeating the checks

From the repository root with the checkout installed:

```bash
python -m pytest -q
python docs/user_manual/render_figures.py
cd docs/user_manual
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

After reviewing the rebuilt PDF, copy it to `PULSUS_User_Manual.pdf` before committing the manual. The publishing steps are recorded in [the release checklist](release_checklist.md).
