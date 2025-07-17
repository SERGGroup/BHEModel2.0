import numpy as np
from scipy.special import kn

__corr = np.array([

    [0.6085, 0.2760, 0.2165],
    [1.3465, 0.4151, 0.7052],
    [0.3777, 0.2792, 0.2195],
    [0.3324, 2.8015, 2.6310],
    [0.2082, 0.0962, 0.5448]

])

def evaluate_ground_f(td, pe=0., use_old_correlation=False):

    if use_old_correlation:
        """
            This function returns the heat transfer between the annulus and the ground [kW/m]
            The implemented model is derived from:

                Y. Zhang, L. Pan, K. Pruess, S. Finsterle, Geothermics (2011), "A time-convolution approach for modeling heat
                exchange between a wellbore and surrounding formation"
                doi: 10.1016/j.geothermics.2011.08.003

        """
        if td < 2.8:

            f = ((np.pi * td) ** -0.5 + 0.5 - (td / np.pi) ** 0.5 / 4 + td / 8)

        else:

            gamma = 0.57722
            theta = (np.log(4 * td) - 2 * gamma)
            f = 2 * (1 / theta - gamma / theta ** 2)

    else:
        """
            This function implements a correlation from FEM results as evaluated 
            by Edoardo Falchini in his PhD thesis (2023).
            
        """
        params = __corr[:, 1] * np.power(pe, __corr[:, 0]) + __corr[:, 2]

        f = params[0]
        for i in range(int(np.floor((len(params) - 1) / 2))):
            f += kn(i, params[1 + i * 2] * np.power(td, params[2 + i * 2]))

    return f