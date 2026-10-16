# PULSUS documentation

These documents describe the public API and examples for PULSUS 0.2.0.

- Start with the [quick start](quickstart.md): install the package and calculate a spectrum from a complete example.
- Read the [function guide and full reference](api.md) for the common functions, inputs, and return values.
- Use [advanced notes](advanced_notes.md) for optional time-transform and numerical details.
- Build the [LaTeX user manual](user_manual/main.tex) for a tutorial and reference PDF.
- Read the [validation record](VALIDATION.md) for the checks and environment used.
- Run the [first tutorial script](../examples/manual/01_first_spectrum.py) from the repository root.

## Where the documentation lives

```text
docs/
  README.md
  quickstart.md
  api.md
  advanced_notes.md
  user_manual/
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

LaTeX is needed only to build the PDF. It is not needed to install PULSUS or run its Python examples. A notebook environment is also optional.

## Updating examples and documentation

1. For documentation maintenance, install the checkout with `python -m pip install -e .` from its root directory.
2. Run `python docs/user_manual/render_figures.py` from the repository root to regenerate PDF figures. User-facing examples display plots with `plt.show()` and do not save files.
3. Edit the corresponding Markdown and LaTeX explanations. Public function signatures should match `src/pulsus/`.
4. Build and inspect the PDF, including every code block and figure.
5. Review `git diff` before committing documentation and tutorial changes.

For changes to scientific implementations, run the applicable tests as well. Installing the `dev` extra and running the complete test suite are useful for contributors; neither is a prerequisite for a user's first spectrum.
