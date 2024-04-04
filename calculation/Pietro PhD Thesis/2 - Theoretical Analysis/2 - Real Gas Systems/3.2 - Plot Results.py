# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np
import copy


# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "0.support", "res"

)

filename = os.path.join(result_folder, "opt_results.npy")
results = np.load(filename)

filename = os.path.join(result_folder, "depths_mesh.npy")
depths = np.load(filename)

filename = os.path.join(result_folder, "grads_mesh.npy")
grads = np.load(filename)

filename = os.path.join(result_folder, "t_rels.npy")
t_ambs = np.load(filename)

res_shape = results.shape
n_points = res_shape[-1]
n_temp = res_shape[1]

# error_max = np.where(results[5] > 1)
# error_zeros = np.where(results[5] <= 0.2)
#
# for i in range(6):
#     results[i][error_max] = np.nan
#     results[i][error_zeros] = np.nan


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


# %%-------------------------------------   PLOT CONTOURS                       -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 4)

i = 0
k = 1
titles = ["Optimal $P_{mod}$ [-]"]
k_mod = [0, 3]
n_levels = 20

if k < 2:

    base_level = np.linspace(1, 5, n_levels)

    if k == 0:
        titles.append("Optimal Work Extraction [kJ/kg]")

    else:
        titles.append("Optimal Exergy Extraction [kJ/kg]")

else:

    titles.append("Optimal Efficiency [-]")
    base_level = np.linspace(

        np.nanmin(results[k, i]),
        np.nanmax(results[k, i]),
        n_levels

    )

levels = [

    base_level,
    np.linspace(

        np.nanmin(results[k + 3, i]),
        np.nanmax(results[k + 3, i]),
        n_levels

    )

]

norm = plt.Normalize(0, 1)
cmap_copy = copy.copy(plt.cm.viridis_r)

for n in range(2):

    cs = axs[n].contourf(depths, grads, results[k + k_mod[n], i], levels=levels[n], extend='both', cmap=cmap_copy)
    fig.colorbar(cs, ax=axs[n], orientation='vertical')

    axs[n].set_xlabel("Depth [m]")
    axs[n].set_ylabel("Gradient [°C/km]")
    axs[n].set_title(titles[n])
    axs[n].set_xscale("log")
    axs[n].set_xlim((500, 5000))
    axs[n].set_ylim((25, 100))

plt.tight_layout(pad=2)
plt.show()
