"""Regenerate PDF figures for maintainers, without file I/O in user examples.

Run from the repository root: python docs/user_manual/render_figures.py
"""

from pathlib import Path
import runpy
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parents[2]
    examples = root / "examples" / "manual"
    directory = root / "docs" / "user_manual" / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(examples))
    scripts = [
        ("01_first_spectrum.py", "first_nr_spectrum.png"),
        ("02_linear_models.py", "linear_models.png"),
        ("03_third_order_models.py", "third_order_models.png"),
        ("04_selected_points.py", "selected_points.png"),
        ("05_time_domain.py", "time_domain.png"),
    ]
    original_show = plt.show
    try:
        plt.show = lambda *args, **kwargs: None
        for script, image in scripts:
            plt.close("all")
            runpy.run_path(str(examples / script), run_name="__main__")
            figures = plt.get_fignums()
            if len(figures) != 1:
                raise RuntimeError(f"Expected one figure from {script}, got {len(figures)}")
            path = directory / image
            plt.figure(figures[0]).savefig(path, dpi=180)
            plt.close("all")
            print("Rendered:", path)
    finally:
        plt.show = original_show
        sys.path.pop(0)


if __name__ == "__main__":
    main()
