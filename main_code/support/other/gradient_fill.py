from matplotlib.patches import Polygon
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

# Function modified from:
#       https://stackoverflow.com/questions/29321835/is-it-possible-to-get-color-gradients-under-curve-in-matplotlib

def gradient_fill(x, y, fill_color=None, ax=None, fill_dir="top", grad_dir="N", **kwargs):

    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.

    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.

    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.

    fill_dir: direction of the filling (can be "top" or "bottom") switch
        the position of the filling to below or above the function

    grad_dir: modify the direction of the gradient
        (accepted "N", "NE", "EN", "E", "SE", "ES", "S", "SW", "WS", "W", "NW", "WN")

    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """

    if ax is None:
        ax = plt.gca()

    line, = ax.plot(x, y, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = __get_z(grad_dir, alpha, fill_color)

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])

    fill_dir = fill_dir.upper()

    if fill_dir == "TOP":

        xy = np.vstack([[xmax, ymax], [xmin, ymax], xy, [xmax, ymax]])

    else:

        xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])

    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)
    ax.autoscale(True)

    return line, im

def __get_z(grad_dir:str, alpha, fill_color, n_pixel=100):

    grad_dir = grad_dir.upper()

    z = np.empty((n_pixel, n_pixel, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)

    z[:, :, :3] = rgb
    alpha_space = np.linspace(alpha, 0, 2 * n_pixel)[:, None]

    for i in range(n_pixel):

        for j in range(n_pixel):

            if grad_dir == "N":
                z[i, j, -1] = alpha_space[2*(n_pixel - i) - 1]

            elif grad_dir == "NE" or grad_dir == "EN":
                z[i, j, -1] = alpha_space[2 * n_pixel - i - j - 1]

            elif grad_dir == "E":
                z[i, j, -1] = alpha_space[2*(n_pixel - j) - 1]

            elif grad_dir == "SE" or grad_dir == "ES":
                z[i, j, -1] = alpha_space[i - j + n_pixel]

            elif grad_dir == "S":
                z[i, j, -1] = alpha_space[2*i]

            elif grad_dir == "SW" or grad_dir == "WS":
                z[i, j, -1] = alpha_space[i + j]

            elif grad_dir == "W":
                z[i, j, -1] = alpha_space[2*j]

            elif grad_dir == "NW" or grad_dir == "WN":
                z[i, j, -1] = alpha_space[-i + j + n_pixel]

            else:

                # by default is N
                z[i, j, -1] = alpha_space[2*(n_pixel - i) - 1]

    return z
