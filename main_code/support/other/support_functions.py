import numpy as np


def get_np_array(min, max, n_elements=10):

    return np.array(range(n_elements)) / (n_elements - 1) * (max - min) + min


def in2m(value_in):
    """
    convert inches to meter
    :param value_in:
    :return: float
    """
    return float(value_in) * 0.0254
