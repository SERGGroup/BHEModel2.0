import numpy as np
from scipy.special import kn

__corr = np.array([

        [0.544987, 0.329228, 0.170987],
        [1.409380, 0.429548, 0.762447],
        [0.416497, 0.291356, 0.165703],
        [0.395478, 2.191704, 1.766232],
        [0.130155, 0.210717, 0.375546],
        [0.343518, 16004.978905, 4.988577],
        [0.040949, 1.101542, 0.356033]


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