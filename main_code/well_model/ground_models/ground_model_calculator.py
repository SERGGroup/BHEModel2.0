import numpy as np

def evaluate_ground_f(td, use_old_correlation=True):

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
            in the excel file that can be found in:

                \\calculation\\Pietro PhD Thesis\\3 - Model Description\\FreeFEM Calculation
                \\2 - Calculations\\res\\results\\1 - base comparison.xlsx

        """
        log_t = np.log(td)
        f = np.exp(0.017749736 * log_t ** 2 -0.316844708 * log_t)

    return f