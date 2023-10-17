# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.power_plants.base_surface_power_plant import BaseSurfacePlant
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   CALCULATIONS OPTIONS                -------------------------------------> #
n_points = 150
n_temp = 5
t_ambs = np.linspace(10, 50, n_temp)
p_rels = np.logspace(-1, 1, n_points)


# %%-------------------------------------   INIT FLUID                          -------------------------------------> #
fluid = "Carbon Dioxide"
bhe_in = PlantThermoPoint([fluid], [1])
bhe_in.set_variable("T", bhe_in.RPHandler.TC)
bhe_in.set_variable("P", bhe_in.RPHandler.PC)
rho_crit = bhe_in.get_variable("rho")


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
depth = 3000        # m
grad_t_rocks = 25   # °C/km

results = np.empty((8, n_temp, n_points))
results[:] = np.nan
elements_list = [list(), list()]

pbar = tqdm(desc="Calculating Points", total=n_points*n_temp)
for i in range(n_temp):

    bhe_in.set_variable("T", t_ambs[i])
    bhe_in.set_variable("rho", rho_crit)
    p_ref = bhe_in.get_variable("P")
    t_rocks = t_ambs[i] + grad_t_rocks * depth / 1000

    elements_list[0].append(list())
    elements_list[1].append(list())

    for j in range(n_points):

        p_curr = p_ref * p_rels[j]
        bhe_in.set_variable("T", t_ambs[i])
        bhe_in.set_variable("P", p_curr)

        bhe = SimplifiedBHE(bhe_in, dz_well=depth, t_rocks=t_rocks)
        base_surface = BaseSurfacePlant(bhe)
        base_surface.prepare_calculation()
        base_surface.calculate()

        elements_list[0][-1].append(bhe)
        elements_list[1][-1].append(base_surface)

        results[0, i, j] = base_surface.w_dot
        results[1, i, j] = base_surface.ex_dot
        results[2, i, j] = base_surface.eta_exs

        results[3, i, j] = base_surface.w_dot_nd
        results[4, i, j] = base_surface.ex_dot_nd
        results[5, i, j] = base_surface.eta_exs

        results[6, i, j] = base_surface.w_dex_min
        results[7, i, j] = base_surface.w_dex_max

        pbar.update(1)

pbar.close()


# %%-------------------------------------   PLOT LINES                          -------------------------------------> #
plot_nd=False

fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

for ax in base_axs[:, 1]:

    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]

if plot_nd:

    k_mod = 3

else:

    k_mod = 0

for i in range(n_temp):

    color = cmap(norm(i / (n_temp - 1)))

    for k in range(len(axs)):
        axs[k].plot(p_rels, results[k + k_mod, i, :], color=color)
        axs[k].set_xscale("log")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   PLOT EXPANSION %                    -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(8, 5)

for i in range(n_temp):

    color = cmap(norm(i / (n_temp - 1)))
    ax.fill_between(p_rels, results[-1, i, :], results[-2, i, :], color=color , alpha=0.1)
    ax.plot(p_rels, results[-1, i, :], color=color, linewidth=2)
    ax.plot(p_rels, results[-2, i, :], color=color, linewidth=2)

ax.set_xscale("log")
plt.tight_layout(pad=2)
plt.show()