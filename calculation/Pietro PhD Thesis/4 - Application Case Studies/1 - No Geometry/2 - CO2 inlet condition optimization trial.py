# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.power_plants.base_surface_power_plant import BaseSurfacePlant
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

def update_values(curr_bhe, t_amb, p_in):

    curr_bhe.input_point.set_variable("T", t_amb)
    curr_bhe.input_point.set_variable("P", p_in)

    new_base_surface = BaseSurfacePlant(curr_bhe)
    new_base_surface.prepare_calculation()
    new_base_surface.calculate()

    return new_base_surface


# %%-------------------------------------   CALCULATIONS OPTIONS                -------------------------------------> #
n_points = 150
n_temp = 5
p_log10_max = 1
t_ambs = np.linspace(10, 50, n_temp)
p_rels = np.logspace(-1, p_log10_max, n_points)


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

pbar = tqdm(desc="Calculating Points", total=n_points*n_temp)
for i in range(n_temp):

    bhe_in.set_variable("T", t_ambs[i])
    bhe_in.set_variable("rho", rho_crit)
    p_ref = bhe_in.get_variable("P")
    t_rocks = t_ambs[i] + grad_t_rocks * depth / 1000
    bhe = SimplifiedBHE(bhe_in, dz_well=depth, t_rocks=t_rocks)

    for j in range(n_points):

        p_curr = p_ref * p_rels[j]
        base_surface = update_values(curr_bhe=bhe, t_amb=t_ambs[i], p_in=p_curr)

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

# %%-------------------------------------   OPTIMIZATION                        -------------------------------------> #
depth = 3000        # m
grad_t_rocks = 25   # °C/km

opt_results = np.empty((2 * 3, n_temp))
opt_results[:] = np.nan

pbar = tqdm(desc="Calculating Points", total=3*n_temp)
for i in range(n_temp):

    bhe_in.set_variable("T", t_ambs[i])
    bhe_in.set_variable("rho", rho_crit)
    p_ref = bhe_in.get_variable("P")
    t_rocks = t_ambs[i] + grad_t_rocks * depth / 1000
    bhe = SimplifiedBHE(bhe_in, dz_well=depth, t_rocks=t_rocks)

    def opt_func(x, opt_value=(0,)):

        p_curr = p_ref * x[0]
        base_surface = update_values(curr_bhe=bhe, t_amb=t_ambs[i], p_in=p_curr)

        if opt_value[0] == 0:

            return -base_surface.w_dot

        elif opt_value[0] == 1:

            return -base_surface.ex_dot

        else:

            return -base_surface.eta_exs

    x0 = [2, 2, 0.8]
    for j in range(3):

        res = minimize(opt_func, np.array(x0[j]), args=[j])
        opt_results[j, i] = res.x[0]
        opt_results[j + 3, i] = -res.fun

        res_test = -opt_func([0.9999], opt_value=(j,))
        if res_test > opt_results[j + 3, i]:
            opt_results[j, i] = 0.9999
            opt_results[j + 3, i] = res_test

        res_test = -opt_func([1.00001], opt_value=(j,))
        if res_test > opt_results[j + 3, i]:
            opt_results[j, i] = 1.00001
            opt_results[j + 3, i] = res_test

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

        if np.log10(opt_results[k, i]) < p_log10_max:

            axs[k].scatter(

                opt_results[k, i], opt_results[k + 3, i], s=30,
                edgecolors=color, facecolors="white", zorder=10

            )

        axs[k].plot(p_rels, results[k + k_mod, i, :], color=color, zorder=1)
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