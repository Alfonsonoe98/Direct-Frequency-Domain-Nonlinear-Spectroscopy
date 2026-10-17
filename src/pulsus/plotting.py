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
    savepath=None,
    dpi=300,
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
    savepath : str or path-like, optional
        Save the figure after plotting. The filename extension selects the
        format, for example PNG or PDF. Default is None, which saves nothing.
        The destination directory must already exist.
    dpi : float, optional
        Resolution used for an optional image export. Default is 300.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the plot.
    ax : matplotlib.axes.Axes
        Axes containing the plot.
    contour : matplotlib.contour.QuadContourSet
        Contour object.

    Notes
    -----
    This function does not call plt.show() or close the figure. Saving exports
    the figure, not the complex numerical response. For further customization,
    modify the returned figure or axes and call fig.savefig(...) afterward.
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

    if savepath is not None:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")

    return fig, ax, contour


def plot_linear_spectrum(
    omega,
    spectrum,
    component="real",
    normalize=False,
    ax=None,
    title=None,
    xlabel=r"$\omega$",
    ylabel=None,
    label=None,
    savepath=None,
    dpi=300,
    **plot_kwargs,
):
    """Plot one component of a calculated complex linear spectrum.

    Parameters
    ----------
    omega : array_like
        Nonempty one-dimensional real frequency axis. The supplied order
        and signs are preserved.
    spectrum : array_like
        One-dimensional complex response, with the same length as omega.
    component : {"real", "imag", "abs"}, optional
        Component to display. Default is "real", matching plot_spectrum.
    normalize : bool, optional
        Divide the displayed component by its maximum absolute value.
        A zero response remains zero. Default is False. The input spectrum
        is never modified.
    ax : matplotlib.axes.Axes, optional
        Existing axes, for example to overlay models. Otherwise create axes.
    title : str, optional
        Plot title. Leave an existing title alone when omitted.
    xlabel, ylabel : str, optional
        Axis labels. By default the vertical label describes the selected
        component and whether it is normalized.
    label : str, optional
        Curve label for a caller-supplied legend.
    savepath : str or path-like, optional
        Save the figure to this path after plotting. The filename extension
        selects the format, for example PNG or PDF. Default is None, which
        saves nothing. The destination directory must already exist.
    dpi : float, optional
        Resolution used for an optional image export. Default is 300.
    **plot_kwargs
        Matplotlib line options, such as color, linestyle, linewidth,
        marker, and markevery. No markers are added by default.

    Returns
    -------
    fig, ax, line
        Matplotlib figure, axes, and Line2D artist, following the
        (figure, axes, artist) convention of pulsus.plot_spectrum.

    Notes
    -----
    This function plots supplied data; it does not calculate a response,
    choose a spectroscopy sign convention, or call plt.show().
    Files are saved only when savepath is supplied. To export a figure
    after adding more curves, legends, or other customization, call
    fig.savefig(...) yourself at the end. Saving exports the figure,
    not the complex numerical response.
    """
    frequency = np.asarray(omega)
    response = np.asarray(spectrum, dtype=complex)
    if frequency.ndim != 1:
        raise ValueError("omega must be one-dimensional")
    if response.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    if frequency.size == 0:
        raise ValueError("omega and spectrum must not be empty")
    if response.shape != frequency.shape:
        raise ValueError(
            "spectrum must have the same shape as omega: "
            f"expected {frequency.shape}, got {response.shape}"
        )
    if np.iscomplexobj(frequency):
        if np.any(frequency.imag != 0):
            raise ValueError("omega must contain real frequencies")
        frequency = frequency.real
    frequency = np.asarray(frequency, dtype=float)
    if not np.isfinite(frequency).all() or not np.isfinite(response).all():
        raise ValueError("omega and spectrum must contain only finite values")
    if not isinstance(component, str):
        raise ValueError("component must be 'real', 'imag', or 'abs'")
    component = component.lower()
    if component == "real":
        data, default_label = response.real, "Real response"
    elif component == "imag":
        data, default_label = response.imag, "Imaginary response"
    elif component == "abs":
        data, default_label = np.abs(response), r"Magnitude $|S|$"
    else:
        raise ValueError("component must be 'real', 'imag', or 'abs'")
    if normalize:
        scale = np.max(np.abs(data))
        if scale > 0:
            data = data / scale
        default_label = (
            r"Normalized magnitude $|S|$"
            if component == "abs"
            else "Normalized " + default_label.lower()
        )
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    plot_kwargs.setdefault("marker", "None")
    line, = ax.plot(frequency, data, label=label, **plot_kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(default_label if ylabel is None else ylabel)
    if title is not None:
        ax.set_title(title)
    if savepath is not None:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")
    return fig, ax, line
