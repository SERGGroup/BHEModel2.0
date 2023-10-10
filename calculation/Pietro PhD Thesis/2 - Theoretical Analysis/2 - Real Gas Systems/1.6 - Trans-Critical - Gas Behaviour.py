# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 500
n_depth = 2
n_v_rel = 10

grad_nd_list = np.logspace(0, 2, n_grad)
dz_nd_list = np.logspace(-2, -1, n_depth)
v_rel_list = np.logspace(0, 2, n_v_rel + 1)[1:]

grad_nd, dz_nd = np.meshgrid(grad_nd_list, dz_nd_list, indexing="ij")

carnot_factor = 1 - 1 / (1 + grad_nd * dz_nd)

spc_work_gas = (grad_nd - 1) * dz_nd
spc_ex_gas = spc_work_gas - np.log((1 + grad_nd * dz_nd)/(1 + dz_nd))
ex_eta_gas = spc_ex_gas / (spc_work_gas * carnot_factor)

spc_work_liq = grad_nd * dz_nd
spc_ex_liq = spc_work_liq - np.log(1 + spc_work_liq)
ex_eta_liq = spc_ex_liq / (spc_work_liq * carnot_factor)

if np.isnan(ex_eta_gas[0, 0]):

    ex_eta_gas[0, :] = np.ones(ex_eta_gas[0, :].shape)

base_shape = np.array(grad_nd.shape)
res_shape = np.append([len(v_rel_list)], base_shape)

r_cp_in_arr = np.empty(len(v_rel_list))
r_cp_in_arr[:] = np.nan

grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
t_rocks = np.empty(res_shape)
w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan
w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
t_rel = 1.5
in_state = PlantThermoPoint(["Methane"], [1], unit_system="MASS BASE SI")
bhe_in = PlantThermoPoint(["Methane"], [1])

in_state.set_variable("P", in_state.RPHandler.PC)
in_state.set_variable("T", in_state.RPHandler.TC)

v_crit = 1 / in_state.get_variable("rho")

states = list()
for i in range(4):
    states.append(in_state.duplicate())

bhe_list = list()
pbar = tqdm(desc="Calculating Points", total=len(v_rel_list) * len(grad_nd_list) * len(dz_nd_list))
for i in range(len(v_rel_list)):

    v_in = v_rel_list[i] * v_crit
    t_in = t_rel * in_state.RPHandler.TC
    in_state.set_variable("T", t_in)
    in_state.set_variable("rho", 1 / v_in)
    in_state.copy_state_to(bhe_in)

    r_cp_in_arr[i] = 1 - 1 / in_state.get_variable("CP/CV")
    grad_rocks[i, :, :] = grad_nd * scipy.constants.g / in_state.get_variable("cp")
    depth[i, :, :] = dz_nd * in_state.get_variable("cp") * t_in / scipy.constants.g
    t_rocks[i, :, :] = t_in + depth[i, :, :] * grad_rocks[i, :, :]

    bhe_sub_list = list()
    for j in range(len(grad_nd_list)):

        bhe_subsub_list = list()
        for k in range(len(dz_nd_list)):

            bhe = SimplifiedBHE(bhe_in, dz_well=depth[i, j, k], t_rocks=t_rocks[i, j, k]-273.15)
            bhe.update()
            bhe_subsub_list.append(bhe)

            for n in range(4):
                bhe.points[n].copy_state_to(states[n])

            surface_states = list()
            surface_states.append(states[3].duplicate())
            surface_states.append(states[3].duplicate())

            surface_states[0].set_variable("P", states[3].get_variable("P"))
            surface_states[0].set_variable("S", states[0].get_variable("S"))

            surface_states[1].set_variable("P", states[0].get_variable("P"))
            surface_states[1].set_variable("S", states[3].get_variable("S"))

            w_dot = states[3].get_variable("H") - states[0].get_variable("H")
            ex_dot = w_dot - states[0].get_variable("T") * (states[3].get_variable("S") - states[0].get_variable("S"))
            ex_dot_nd = ex_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))

            if (not (w_dot < 0 or ex_dot < 0)) and (not states[-1].get_variable("H") < -1e6) and (grad_nd[j, k] < 25):

                cf = 1 - states[0].get_variable("T") / states[2].get_variable("T")

                w_dot_nds[i, j, k] = w_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                ex_dot_nds[i, j, k] = ex_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                eta_exs[i, j, k] = ex_dot / (w_dot * cf)

                w_dex_mins[i, j, k] = (surface_states[1].get_variable("H") - states[0].get_variable("H")) / w_dot
                w_dex_maxs[i, j, k] = (states[3].get_variable("H") - surface_states[1].get_variable("H")) / w_dot

            pbar.update(1)

        bhe_sub_list.append(bhe_subsub_list)

    bhe_list.append(bhe_sub_list)

pbar.close()


# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

for ax in base_axs[:, 1]:

    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]
x_values = grad_nd - 1

v_rel_label = "$v_{{rel}} = 10^{{ {:0.1f} }}$"
cmp_labels = ["Liquid", v_rel_label, "Ideal Gas"]
cmp_y_values = [

    [spc_work_liq, spc_ex_liq, ex_eta_liq],
    [w_dot_nds, ex_dot_nds, eta_exs],
    [spc_work_gas, spc_ex_gas, ex_eta_gas],

]

k = 0
lines = [list(), list(), list()]

for m in range(len(cmp_y_values[0])):

    for n in range(len(cmp_labels)):

        if cmp_labels[n] == v_rel_label:

            for i in range(len(v_rel_list)):

                lines[m].append(

                    axs[m].plot(

                        x_values[:, k], cmp_y_values[n][m][i, :, k], "-",
                        label=v_rel_label.format(np.log10(v_rel_list[i])),
                        color=cmap(norm((i + 1) / (len(v_rel_list) + 1)))

                    )[0]

                )

        else:

            lines[m].append(

                axs[m].plot(

                    x_values[:, k], cmp_y_values[n][m][:, k], "-",
                    label=cmp_labels[n], color=cmap(norm(n)),
                    linewidth=2

                )[0]

            )


y_names = ["${\\dot{w}}^{\\#}$ [-]", "${\\dot{e}_x}^{\\#}$ [-]", "${\\eta}_{ex}$ [-]"]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#} - 1$ [-]")

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")

    if not y_names[k] == y_names[1]:
        axs[k].legend(handles=lines[k], fontsize="8")

axs[-1].set_ylim((

    0.3,
    1

))

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   EVALUATE COMPRESSIBILITY EFFECT     -------------------------------------> #
t_in = t_rel * in_state.RPHandler.TC
p_in = 0.01 * in_state.RPHandler.PC
in_state.set_variable("T", t_in)
in_state.set_variable("P", p_in)
ig_comp = 1 - 1 / in_state.get_variable("CP/CV")

k = 0
min_grad_nd = list()
min_grad_calc = list()
comp_rel = r_cp_in_arr

for i in range(len(v_rel_list)):

    arr = w_dot_nds[i, :, k]
    min_index = np.where(arr == np.nanmin(arr))[0]
    min_grad_nd.append(grad_nd[min_index[0], k])

    v_in = v_rel_list[i] * v_crit
    t_in = t_rel * in_state.RPHandler.TC
    in_state.set_variable("T", t_in)
    in_state.set_variable("rho", 1 / v_in)

    dpdt_in = in_state.get_derivative("P", "T", "rho")
    gamma_in = in_state.get_variable("CP/CV")
    cp_in = in_state.get_variable("cp")
    alpha_in = cp_in * (1 - 1 / gamma_in) / (v_in * dpdt_in)
    min_grad_calc.append(alpha_in)

n = 90
plt.plot(comp_rel, min_grad_nd)
plt.scatter(comp_rel, min_grad_calc)
plt.scatter(ig_comp, 1)
plt.xscale("log")
# plt.yscale("log")
plt.show()
