import numpy as np


def get_np_array(range_min, range_max, n_elements=10, log_scale=False):

    if log_scale:
        return np.power(range_max / range_min, np.array(range(n_elements)) / (n_elements - 1)) * range_min

    else:
        return np.array(range(n_elements)) / (n_elements - 1) * (range_max - range_min) + range_min


def in2m(value_in):
    """
    convert inches to meter
    :param value_in:
    :return: float
    """
    return float(value_in) * 0.0254
