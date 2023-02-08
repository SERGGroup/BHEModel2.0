import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# Function modified from: https://matplotlib.org/stable/gallery/lines_bars_and_markers/curve_error_band.html
# using: https://www.geeksforgeeks.org/matplotlib-pyplot-fill_between-in-python/
def draw_error_band(ax, x, y, err, err_fixed=False, alpha=0., **kwargs):

    if err_fixed:

        y_plus = y + err
        y_minus = y - err

    else:

        y_plus = y*(1 + err)
        y_minus = y * (1 - err)

    line = ax.plot(x, y,  **kwargs)
    ax.fill_between(x, y_plus, y_minus, alpha=alpha, color=line[0].get_color(), **kwargs)

    return line[0]