# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 25
n_depth = 3
n_t_rel = 5

grad_nd_list = np.logspace(0, 2, n_grad + 1)[1:]
dz_nd_list = np.logspace(-2, 0, n_depth)
t_rel_list = np.logspace(np.log10(0.5), np.log10(8), n_t_rel)
grad_nd, dz_nd = np.meshgrid(grad_nd_list, dz_nd_list, indexing="ij")

base_shape = np.array(grad_nd.shape)
res_shape = np.append([len(t_rel_list)], base_shape)

grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
t_rocks = np.empty(res_shape)

carnot_factor = 1 - 1 / (1 + grad_nd * dz_nd)

spc_work_gas = (grad_nd - 1) * dz_nd
spc_ex_gas = spc_work_gas - np.log((1 + grad_nd * dz_nd)/(1 + dz_nd))
ex_eta_gas = spc_ex_gas / (spc_work_gas * carnot_factor)

spc_work_liq = grad_nd * dz_nd
spc_ex_liq = spc_work_liq - np.log(1 + spc_work_liq)
ex_eta_liq = spc_ex_liq / (spc_work_liq * carnot_factor)

grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan

w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
p_rel_curr = 10**3
in_state = PlantThermoPoint(["Methane"], [1], unit_system="MASS BASE SI")
bhe_in = PlantThermoPoint(["Methane"], [1])

states = list()
for i in range(4):
    states.append(in_state.duplicate())

bhe_list = list()
p_in = in_state.RPHandler.PC * p_rel_curr
pbar = tqdm(desc="Calculating Points", total=len(t_rel_list) * len(grad_nd_list) * len(dz_nd_list))
for i in range(len(t_rel_list)):

    t_rel = t_rel_list[i]
    t_in = t_rel * in_state.RPHandler.TC
    in_state.set_variable("T", t_in)
    in_state.set_variable("P", p_in)
    in_state.copy_state_to(bhe_in)

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

            if not (w_dot < 0 or ex_dot < 0):
                w_dot_nds[i, j, k] = w_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                ex_dot_nds[i, j, k] = ex_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                eta_exs[i, j, k] = ex_dot / (w_dot * carnot_factor[j, k])

                w_dex_mins[i, j, k] = (surface_states[0].get_variable("H") - states[0].get_variable("H")) / w_dot
                w_dex_maxs[i, j, k] = (states[3].get_variable("H") - surface_states[1].get_variable("H")) / w_dot

            pbar.update(1)

        bhe_sub_list.append(bhe_subsub_list)

    bhe_list.append(bhe_sub_list)

pbar.close()


# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
colors = [

    "tab:blue", "tab:orange",
    "tab:green", "tab:red",
    "tab:purple"

]

list_styles = ["--", "-.", "-", "x", "."]

label_str = "${{\\Delta z}}^{{\\#}}\\ =\\ 10^{{ {order} }}$"
fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
# remove the underlying axes
for ax in base_axs[:, 1]:
    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]

x_values = grad_nd - 1

for i in range(len(t_rel_list)):

    lines = [list(), list(), list()]

    for k in range(len(dz_nd_list)):

        color = colors[k]
        dz_nd_curr = dz_nd_list[k]
        label_curr = label_str.format(

            order=np.round(np.log10(dz_nd_curr), 1)

        )

        lines[0].append(

            axs[0].plot(

                x_values[:, k],
                w_dot_nds[i, :, k],
                list_styles[i], label=label_curr,
                color=color

            )[0]

        )

        lines[1].append(

            axs[1].plot(

                x_values[:, k],
                ex_dot_nds[i, :, k],
                list_styles[i], label=label_curr,
                color=color

            )[0]

        )

        lines[2].append(

            axs[2].plot(

                x_values[:, k],
                eta_exs[i, :, k],
                list_styles[i], label=label_curr,
                color=color

            )[0]

        )

y_names = ["${\\dot{w}}^{\\#}$ [-]", "${\\dot{e}_x}^{\\#}$ [-]", "${\\eta}_{ex}$ [-]"]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#} - 1$ [-]")

    # axs[k].set_xlim((
    #
    #     np.min(x_values),
    #     np.max(x_values)
    #
    # ))

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")
        # axs[k].set_ylim((
        #
        #     np.nanmin(ex_dot_nds) / 2,
        #     np.nanmin(ex_dot_nds) * 2
        #
        # ))

    if not y_names[k] == y_names[1]:
        axs[k].legend(handles=lines[k], fontsize="8")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()
