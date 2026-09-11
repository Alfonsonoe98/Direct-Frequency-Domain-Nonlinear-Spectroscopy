import numpy as np
import matplotlib.pyplot as plt


def plot_spectrum(
    omega1,
    omega3,
    spectrum,
    component="real",
    normalize=False,
    levels=31,
    ax=None,
    colorbar=True,
    colorbar_label=None,
    title=None,
    xlabel=r"$\omega_1$",
    ylabel=r"$\omega_3$",
):
    """
    Plot a two-dimensional complex spectrum.

    Parameters
    ----------
    omega1 : array_like
        Excitation-frequency axis.
    omega3 : array_like
        Detection-frequency axis.
    spectrum : array_like
        Complex spectrum with shape
        (len(omega3), len(omega1)).
    component : {"real", "imag", "abs"}, optional
        Component of the complex spectrum to display.
    normalize : bool, optional
        Normalize the displayed component by its maximum
        absolute value. Default is False.
    levels : int, optional
        Number of contour levels.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to plot.
    colorbar : bool, optional
        Add a colorbar.
    colorbar_label : str, optional
        Custom colorbar label. If omitted, PULSUS chooses
        a label from component and normalize.
    title : str, optional
        Plot title.
    xlabel, ylabel : str, optional
        Axis labels.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the plot.
    ax : matplotlib.axes.Axes
        Axes containing the plot.
    contour : matplotlib.contour.QuadContourSet
        Contour object.
    """
    omega1 = np.asarray(
        omega1,
        dtype=float,
    )

    omega3 = np.asarray(
        omega3,
        dtype=float,
    )

    spectrum = np.asarray(
        spectrum,
        dtype=complex,
    )

    if omega1.ndim != 1 or omega3.ndim != 1:
        raise ValueError(
            "omega1 and omega3 must be one-dimensional"
        )

    expected_shape = (
        omega3.size,
        omega1.size,
    )

    if spectrum.shape != expected_shape:
        raise ValueError(
            "spectrum must have shape "
            f"{expected_shape}, got {spectrum.shape}"
        )

    component = component.lower()

    if component == "real":
        data = spectrum.real
        default_label = "Real response"

    elif component == "imag":
        data = spectrum.imag
        default_label = "Imaginary response"

    elif component == "abs":
        data = np.abs(spectrum)
        default_label = r"Magnitude $|S|$"

    else:
        raise ValueError(
            "component must be 'real', 'imag', or 'abs'"
        )

    if normalize:
        scale = np.max(
            np.abs(data)
        )

        if scale > 0:
            data = data / scale

        if component == "abs":
            default_label = r"Normalized magnitude $|S|$"
        else:
            default_label = (
                "Normalized "
                + default_label.lower()
            )

    if colorbar_label is None:
        colorbar_label = default_label

    if ax is None:
        fig, ax = plt.subplots()

    else:
        fig = ax.figure

    if component in ("real", "imag"):
        vmax = np.max(
            np.abs(data)
        )

        if vmax == 0:
            vmax = 1.0

        contour_levels = np.linspace(
            -vmax,
            vmax,
            levels,
        )

        contour = ax.contourf(
            omega1,
            omega3,
            data,
            levels=contour_levels,
            cmap="RdBu_r",
        )

    else:
        vmax = np.max(data)

        if vmax == 0:
            vmax = 1.0

        contour_levels = np.linspace(
            0.0,
            vmax,
            levels,
        )

        contour = ax.contourf(
            omega1,
            omega3,
            data,
            levels=contour_levels,
            cmap="viridis",
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    if colorbar:
        cbar = fig.colorbar(
            contour,
            ax=ax,
        )

        cbar.set_label(
            colorbar_label
        )

    return fig, ax, contour
