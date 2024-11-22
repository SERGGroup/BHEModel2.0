# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.ndimage import gaussian_filter, uniform_filter
from matplotlib.animation import FuncAnimation
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import pandas as pd
import numpy as np
import warnings
import os


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
output_folder = os.path.join(base_folder, "00 - Output", "3 - Overall Optimization Results", "2024-10-29 - Main Results")

opt_res = np.load(os.path.join(output_folder, 'optimization_results.npy'))
x, y = np.meshgrid(opt_res[:, 0, 0, 0], opt_res[0, :, 0, 1], indexing="ij")


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
omega_i = 80
plot_mean = True
window_size = 9
omega_range = np.intc(np.linspace(-window_size, window_size, 2*window_size+1))

w_net = opt_res[:, :, 0, 5]
m_dot = opt_res[:, :, 0, 6]
p_max = opt_res[:, :, 0, 7]
t_ihx = opt_res[:, :, 0, 8]
ihx_power = opt_res[:, :, 0, 9]

omega = opt_res[0, 0, omega_i, 2]

if plot_mean:
    omega_i += omega_range

for i in range(len(x)):

    for j in range(len(y)):

        w_net[i, j] = np.nanmean(opt_res[i, j, omega_i, 5])
        m_dot[i, j] = np.nanmean(opt_res[i, j, omega_i, 6])
        p_max[i, j] = np.nanmean(opt_res[i, j, omega_i, 7])
        t_ihx[i, j] = np.nanmean(opt_res[i, j, omega_i, 8])
        ihx_power[i, j] = np.nanmean(opt_res[i, j, omega_i, 9])

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# <-- FIRST AX ----------------------------------->
axs[0].set_title("Power Requirement [-]")

contour = axs[0].contourf(y, x, w_net)
fig.colorbar(contour, ax=axs[0])

# <-- SECOND AX ---------------------------------->
axs[1].set_title("m dot rel [-]")
contour = axs[1].contourf(y, x, m_dot)
fig.colorbar(contour, ax=axs[1])

for ax in axs:

    ax.set_xlabel("Geothermal Gradient [°C/km]")
    ax.set_ylabel("Depth [km]")
    ax.invert_yaxis()

fig.suptitle("Omega = {:.2f}".format(omega))
plt.tight_layout(pad=1)
plt.show()
