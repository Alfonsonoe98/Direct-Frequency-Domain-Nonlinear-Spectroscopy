# PULSUS documentation

These documents describe **PULSUS 0.3.0**. The linear plotter and optional figure-export arguments require 0.3.0 or newer; an older installation must be upgraded or replaced by this checkout. The first quick-start calculation also works with 0.2.0.

- Start with the [quick start](quickstart.md): install the package and calculate a spectrum from a complete example.
- Read the [function guide and full reference](api.md) for the common functions, inputs, and return values.
- Use [method and scientific validation notes](method_notes.md) for the pulse-dressing formulation, conventions, and benchmark results.
- Use [advanced notes](advanced_notes.md) for optional time-transform and numerical details.
- Read the [user manual (PDF)](user_manual/PULSUS_User_Manual.pdf) for tutorials and reference material. The [LaTeX source](user_manual/main.tex) is available for documentation maintenance.
- Read the [documentation and package validation record](VALIDATION.md) for software checks and the environment used.
- Follow the [release checklist](release_checklist.md) to merge and publish version 0.3.0.
- Run the [first tutorial script](../examples/manual/01_first_spectrum.py) from the repository root.

## Where the documentation lives

```text
docs/
  README.md
  quickstart.md
  api.md
  advanced_notes.md
  method_notes.md
  release_checklist.md
  releases/
    0.3.0.md              # GitHub release notes
  user_manual/
    PULSUS_User_Manual.pdf  # published, ready-to-read manual
    main.tex
    render_figures.py      # maintainer-only figure generation
    sections/
    figures/
    build/                 # generated files; do not commit
examples/
  manual/
    01_first_spectrum.py   # copy/pasteable; displays a plot, saves no files
    02_linear_models.py
    03_third_order_models.py
    04_selected_points.py
    05_time_domain.py
    common.py              # shared model for examples 02--05
```

The GitHub repository is the canonical home of the manual. Overleaf is optional: local editing and compilation are sufficient. Keep tutorial programs under `examples/manual/`, run them when updating the documentation, and include them in LaTeX with `\lstinputlisting` so the PDF uses the same code as the runnable example.

## Build the PDF locally

Install a TeX distribution with `latexmk` and `pdflatex`. These are maintainer steps; ordinary PULSUS users need only the pip package. From the repository root:

```bash
python docs/user_manual/render_figures.py
cd docs/user_manual
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The resulting PDF is `docs/user_manual/build/main.pdf`. On macOS, open it with `open build/main.pdf`; on Linux use `xdg-open build/main.pdf`; on Windows open the file in your PDF viewer. The `build/` directory is ignored by Git.

After reviewing the rebuilt PDF, refresh the published copy from `docs/user_manual/`:

```bash
cp build/main.pdf PULSUS_User_Manual.pdf
```

Commit this published PDF with any changes to the manual source or figures. Keep the temporary files in `build/` untracked.

LaTeX is needed only to build the PDF. It is not needed to install PULSUS or run its Python examples. A notebook environment is also optional.

## Updating examples and documentation

1. For documentation maintenance, install the checkout with `python -m pip install -e .` from its root directory so examples use the same code as the documentation.
2. Run `python docs/user_manual/render_figures.py` from the repository root to regenerate PDF figures. User-facing examples display plots with `plt.show()` and do not save files.
3. Edit the corresponding Markdown and LaTeX explanations. Public function signatures should match `src/pulsus/`.
4. Build and inspect the PDF, including every code block and figure.
5. Copy the reviewed `build/main.pdf` to `PULSUS_User_Manual.pdf`.
6. Review `git diff` before committing the documentation, tutorials, and published PDF together.

For changes to scientific implementations, run the applicable tests as well. Installing the `dev` extra and running the complete test suite are useful for contributors; neither is a prerequisite for a user's first spectrum.

Both plotting helpers return Matplotlib figures and axes. User examples display them with `plt.show()` and save nothing by default. Optional `savepath` and `dpi` arguments export a PNG or PDF figure; `fig.savefig(...)` remains available for saving after further modifications. Keep version information in the README, quick start, reference, and PDF synchronized when preparing a release.
