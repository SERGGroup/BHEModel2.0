# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "4 - HTHP",
    "res"

)

filename = os.path.join(result_folder, "opt_results.npy")
results = np.load(filename)

filename = os.path.join(result_folder, "depths_mesh.npy")
depths = np.load(filename)

filename = os.path.join(result_folder, "grads_mesh.npy")
grads = np.load(filename)

filename = os.path.join(result_folder, "t_ambs.npy")
t_ambs = np.load(filename)

res_shape = results.shape
n_points = res_shape[-1]
n_temp = res_shape[1]


# %%-------------------------------------   PLOT LINES                          -------------------------------------> #
fig, axs = plt.subplots(1, 3, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

i = 1
for j in range(n_points):

    color = cmap(norm(j / (n_points - 1)))

    for k in range(len(axs)):

        axs[k].plot(depths[:, j], results[k, i, :, j], color=color, zorder=1)
        axs[k].set_xscale("log")

        if k < len(axs) - 1:
            axs[k].set_yscale("log")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()

