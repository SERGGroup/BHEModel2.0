# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "1 - No Geometry",
    "res", "support"

)
filename = os.path.join(result_folder, "Optimization Result.npy")
results = np.load(filename)
res_shape = results.shape

n_points = res_shape[-1]
n_temp = res_shape[1]

depths_list = np.logspace(2, 4, n_points)
grads_list = np.logspace(1, 2, n_points)
depths, grads = np.meshgrid(depths_list, grads_list, indexing="ij")

error_max = np.where(results[5] > 1)
error_zeros = np.where(results[5] <= 0.2)

for i in range(6):
    results[i][error_max] = np.nan
    results[i][error_zeros] = np.nan


# %%-------------------------------------   PLOT LINES                          -------------------------------------> #
fig, axs = plt.subplots(1, 3, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

i = 2
for j in range(n_points):

    color = cmap(norm(j / (n_points - 1)))

    for k in range(len(axs)):

        axs[k].plot(depths_list, results[k, i, :, j], color=color, zorder=1)
        axs[k].set_xscale("log")

        if k < len(axs) - 1:
            axs[k].set_yscale("log")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   PLOT CONTOURS                       -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 5)

i = 1
k = 0
k_mod = [0, 3]
for n in range(2):

    axs[n].contourf(depths, grads, results[k + k_mod[n], i], extend='both')
    axs[n].set_xscale("log")
    axs[n].set_xlim((200, 10000))
    axs[n].set_ylim((25, 100))

plt.show()