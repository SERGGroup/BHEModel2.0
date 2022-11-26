import matplotlib.colors as mcolors
import numpy as np

BASE_COLORS = ["#1f77b4", "#ff7f0e"]


def format_title_from_key(key):

    if not "\n" in key:

        return r"$\bf{" + key.replace(" ", "\ ") + "}$" + "\n"

    else:

        spl_key = key.split("\n")
        return r"$\bf{" + spl_key[0].replace(" ", "\ ") + "}$" + "\n" + spl_key[1].replace(" ", "\ ")


class ColorFader:

    def __init__(self, from_value=0, to_value=1):

        self.from_value = from_value
        self.delta_value = to_value - from_value

        self.c0 = np.array(mcolors.to_rgb(BASE_COLORS[0]))
        self.c1 = np.array(mcolors.to_rgb(BASE_COLORS[1]))

    def get_color(self, x):

        perc = (x - self.from_value) / self.delta_value
        return mcolors.to_hex(self.c1 * (1 - perc) + self.c0 * perc)
