import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pulsus


def test_plot_spectrum_automatic_colorbar_label():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.array(
        [
            [1.0 + 1.0j, 0.5 + 0.2j, -0.5 + 0.4j],
            [0.2 - 0.3j, 0.8 + 0.7j, -0.2 + 0.1j],
        ],
        dtype=complex,
    )

    fig, ax, contour = pulsus.plot_spectrum(
        omega1=omega1,
        omega3=omega3,
        spectrum=spectrum,
        component="imag",
        normalize=True,
    )

    colorbar_axis = fig.axes[-1]

    assert colorbar_axis.get_ylabel() == (
        "Normalized imaginary response"
    )

    plt.close(fig)


def test_plot_spectrum_custom_colorbar_label():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.ones(
        (2, 3),
        dtype=complex,
    )

    fig, ax, contour = pulsus.plot_spectrum(
        omega1=omega1,
        omega3=omega3,
        spectrum=spectrum,
        component="real",
        colorbar_label=r"$\mathrm{Re}\,P^{(3)}$",
    )

    colorbar_axis = fig.axes[-1]

    assert colorbar_axis.get_ylabel() == (
        r"$\mathrm{Re}\,P^{(3)}$"
    )

    plt.close(fig)

def test_plot_spectrum():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.array([
        [1.0 + 0.2j, -0.5 + 0.1j, 0.2 - 0.1j],
        [-0.8 + 0.3j, 0.4 - 0.2j, 0.1 + 0.5j],
    ])

    fig, ax, contour = pulsus.plot_spectrum(
        omega1=omega1,
        omega3=omega3,
        spectrum=spectrum,
        component="real",
    )

    assert fig is not None
    assert ax is not None
    assert contour is not None

    assert ax.get_xlabel() == r"$\omega_1$"
    assert ax.get_ylabel() == r"$\omega_3$"

    plt.close(fig)


def test_plot_spectrum_accepts_components():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    spectrum = np.ones(
        (2, 3),
        dtype=complex,
    )

    for component in (
        "real",
        "imag",
        "abs",
    ):
        fig, ax, contour = pulsus.plot_spectrum(
            omega1=omega1,
            omega3=omega3,
            spectrum=spectrum,
            component=component,
            normalize=True,
        )

        assert contour is not None

        plt.close(fig)


def test_plot_spectrum_rejects_wrong_shape():
    omega1 = np.array([
        0.9,
        1.0,
        1.1,
    ])

    omega3 = np.array([
        1.0,
        1.2,
    ])

    wrong_shape = np.zeros(
        (3, 2),
        dtype=complex,
    )

    with pytest.raises(ValueError):
        pulsus.plot_spectrum(
            omega1=omega1,
            omega3=omega3,
            spectrum=wrong_shape,
        )


@pytest.fixture(autouse=True)
def close_plot_figures():
    yield
    plt.close("all")


@pytest.fixture
def linear_data():
    # Preserve negative frequencies and a descending axis; no sorting or
    # spectroscopy sign conversion belongs in a plotting helper.
    return (
        np.array([1.5, 0.2, -0.7]),
        np.array([-2.0 + 5.0j, 1.0 - 2.0j, 0.0 + 1.0j]),
    )


@pytest.fixture
def grid_data():
    return (
        np.array([0.9, 1.0, 1.1]),
        np.array([1.0, 1.2]),
        np.array([
            [1.0 + 0.2j, -0.5 + 0.1j, 0.2 - 0.1j],
            [-0.8 + 0.3j, 0.4 - 0.2j, 0.1 + 0.5j],
        ]),
    )


def call_plotter(kind, linear_data, grid_data, **kwargs):
    if kind == "linear":
        return pulsus.plot_linear_spectrum(*linear_data, **kwargs)
    return pulsus.plot_spectrum(*grid_data, **kwargs)


def test_linear_plotter_is_public():
    from pulsus.plotting import plot_linear_spectrum

    assert pulsus.plot_linear_spectrum is plot_linear_spectrum
    assert "plot_linear_spectrum" in pulsus.__all__


@pytest.mark.parametrize(
    "component, expected, ylabel",
    [
        ("real", [-2.0, 1.0, 0.0], "Real response"),
        ("IMAG", [5.0, -2.0, 1.0], "Imaginary response"),
        ("abs", [np.sqrt(29.0), np.sqrt(5.0), 1.0], r"Magnitude $|S|$"),
    ],
)
def test_linear_components_preserve_supplied_values(
    linear_data, component, expected, ylabel,
):
    fig, ax, line = pulsus.plot_linear_spectrum(
        *linear_data, component=component,
    )
    np.testing.assert_array_equal(line.get_xdata(), linear_data[0])
    np.testing.assert_allclose(line.get_ydata(), expected)
    assert fig is ax.figure
    assert line in ax.lines
    assert ax.get_xlabel() == r"$\omega$"
    assert ax.get_ylabel() == ylabel


def test_linear_defaults_to_real_without_normalization(linear_data):
    _, _, line = pulsus.plot_linear_spectrum(*linear_data)
    np.testing.assert_array_equal(line.get_ydata(), [-2.0, 1.0, 0.0])


@pytest.mark.parametrize(
    "component, expected, ylabel",
    [
        ("real", [-1.0, 0.5, 0.0], "Normalized real response"),
        ("imag", [1.0, -0.4, 0.2], "Normalized imaginary response"),
        (
            "abs",
            [1.0, np.sqrt(5.0 / 29.0), 1.0 / np.sqrt(29.0)],
            r"Normalized magnitude $|S|$",
        ),
    ],
)
def test_linear_normalizes_selected_component_without_changing_input(
    linear_data, component, expected, ylabel,
):
    omega, spectrum = linear_data
    original_omega, original_spectrum = omega.copy(), spectrum.copy()
    omega.flags.writeable = False
    spectrum.flags.writeable = False
    _, ax, line = pulsus.plot_linear_spectrum(
        omega, spectrum, component=component, normalize=True,
    )
    np.testing.assert_allclose(line.get_ydata(), expected)
    np.testing.assert_array_equal(omega, original_omega)
    np.testing.assert_array_equal(spectrum, original_spectrum)
    assert ax.get_ylabel() == ylabel


@pytest.mark.parametrize("component", ["real", "imag", "abs"])
def test_linear_zero_spectrum_remains_finite_when_normalized(component):
    _, _, line = pulsus.plot_linear_spectrum(
        [0.0, 1.0], [0.0j, 0.0j], component=component, normalize=True,
    )
    np.testing.assert_array_equal(line.get_ydata(), [0.0, 0.0])


def test_linear_reuses_axes_and_supports_custom_labels_styles(linear_data):
    fig, ax = plt.subplots()
    ax.set_title("Existing title")
    result_fig, result_ax, first = pulsus.plot_linear_spectrum(
        *linear_data, ax=ax, label="First", color="purple", linestyle="--",
        linewidth=2.5, xlabel="Frequency", ylabel="Response",
    )
    _, _, second = pulsus.plot_linear_spectrum(
        *linear_data, ax=ax, component="imag", label="Second", title="Models",
        xlabel="Frequency", ylabel="Response", marker="o",
    )
    assert result_fig is fig
    assert result_ax is ax
    assert list(ax.lines) == [first, second]
    assert first.get_label() == "First"
    assert first.get_color() == "purple"
    assert first.get_linestyle() == "--"
    assert first.get_linewidth() == 2.5
    assert second.get_marker() == "o"
    assert ax.get_title() == "Models"
    assert ax.get_xlabel() == "Frequency"
    assert ax.get_ylabel() == "Response"


def test_linear_keeps_existing_title_when_title_is_omitted(linear_data):
    fig, ax = plt.subplots()
    ax.set_title("Existing title")
    pulsus.plot_linear_spectrum(*linear_data, ax=ax)
    assert ax.get_title() == "Existing title"


def test_linear_has_no_default_markers_even_with_matplotlib_marker_style(linear_data):
    with matplotlib.rc_context({"lines.marker": "o"}):
        _, _, line = pulsus.plot_linear_spectrum(*linear_data)
    assert line.get_marker() == "None"


@pytest.mark.parametrize(
    "omega, spectrum, message",
    [
        (1.0, [1.0], "omega must be one-dimensional"),
        ([[1.0, 2.0]], [1.0, 2.0], "omega must be one-dimensional"),
        ([1.0], 1.0, "spectrum must be one-dimensional"),
        ([1.0, 2.0], [[1.0, 2.0]], "spectrum must be one-dimensional"),
        ([], [], "must not be empty"),
        ([1.0, 2.0], [1.0], "same shape as omega"),
        ([1.0 + 0.1j], [1.0], "real frequencies"),
        ([np.nan], [1.0], "finite values"),
        ([np.inf], [1.0], "finite values"),
        ([1.0], [np.nan + 0.0j], "finite values"),
        ([1.0], [0.0 + np.inf * 1j], "finite values"),
    ],
)
def test_linear_rejects_invalid_data(omega, spectrum, message):
    with pytest.raises(ValueError, match=message):
        pulsus.plot_linear_spectrum(omega, spectrum)


@pytest.mark.parametrize("component", ["phase", "", None, 1])
def test_linear_rejects_invalid_components(linear_data, component):
    with pytest.raises(ValueError, match="component must be"):
        pulsus.plot_linear_spectrum(*linear_data, component=component)


def test_linear_accepts_real_frequencies_stored_as_complex():
    _, _, line = pulsus.plot_linear_spectrum([2.0 + 0.0j, 1.0], [1.0, 2.0])
    np.testing.assert_array_equal(line.get_xdata(), [2.0, 1.0])


@pytest.mark.parametrize("kind", ["linear", "2d"])
@pytest.mark.parametrize("extension, signature", [("png", b"\x89PNG\r\n\x1a\n"), ("pdf", b"%PDF-")])
def test_plotters_optionally_export_png_and_pdf_and_remain_editable(
    kind, extension, signature, tmp_path, linear_data, grid_data,
):
    fig, ax = plt.subplots()
    destination = tmp_path / f"{kind}.{extension}"
    plotted_fig, plotted_ax, artist = call_plotter(
        kind, linear_data, grid_data, ax=ax, savepath=destination, dpi=120,
    )
    assert plotted_fig is fig
    assert plotted_ax is ax
    assert artist is not None
    assert destination.read_bytes().startswith(signature)
    assert plt.fignum_exists(fig.number)

    # Optional export must preserve the returned Matplotlib objects so callers
    # can continue editing and export their finished figure themselves.
    ax.set_title("Changed after plotting")
    ax.axhline(0.0, color="black")
    final_path = tmp_path / f"{kind}-modified.{extension}"
    fig.savefig(final_path, dpi=90)
    assert final_path.read_bytes().startswith(signature)
    assert ax.get_title() == "Changed after plotting"


@pytest.mark.parametrize("kind", ["linear", "2d"])
def test_plotters_do_not_save_show_or_close_by_default(
    kind, linear_data, grid_data, monkeypatch,
):
    from matplotlib.figure import Figure

    def unexpected_call(*args, **kwargs):
        pytest.fail("Plotting must not save, show, or close by default")

    # Create the axes before patching so Matplotlib backend initialization is
    # outside the function behavior being checked.
    fig, ax = plt.subplots()
    with monkeypatch.context() as patches:
        patches.setattr(Figure, "savefig", unexpected_call)
        patches.setattr(plt, "show", unexpected_call)
        patches.setattr(plt, "close", unexpected_call)
        plotted_fig, _, _ = call_plotter(
            kind, linear_data, grid_data, ax=ax,
        )
    assert plotted_fig is fig
    assert plt.fignum_exists(fig.number)


@pytest.mark.parametrize("kind", ["linear", "2d"])
@pytest.mark.parametrize("dpi", [None, 144])
def test_plotters_export_requested_resolution(
    kind, dpi, linear_data, grid_data, tmp_path,
):
    from PIL import Image

    fig, ax = plt.subplots()
    path = tmp_path / "spectrum.png"
    options = {} if dpi is None else {"dpi": dpi}
    call_plotter(kind, linear_data, grid_data, ax=ax, savepath=path, **options)
    with Image.open(path) as image:
        np.testing.assert_allclose(
            image.info["dpi"], [300 if dpi is None else dpi] * 2, atol=0.05,
        )


def test_2d_plotter_preserves_original_positional_arguments(grid_data):
    fig, ax = plt.subplots()
    actual_fig, actual_ax, contour = pulsus.plot_spectrum(
        *grid_data, "imag", True, 21, ax, False, "Unused", "Legacy call", "X", "Y",
    )
    assert actual_fig is fig
    assert actual_ax is ax
    assert contour is not None
    assert ax.get_title() == "Legacy call"
    assert ax.get_xlabel() == "X"
    assert ax.get_ylabel() == "Y"
    assert len(fig.axes) == 1
