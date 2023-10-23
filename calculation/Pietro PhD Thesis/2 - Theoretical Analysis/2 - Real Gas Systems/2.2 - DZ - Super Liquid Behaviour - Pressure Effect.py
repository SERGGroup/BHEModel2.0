# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import numpy as np
import gc


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 350
n_depth = 5

grad_nd_nd_list = np.logspace(-2, 4.5, n_grad)
dz_nd_list = np.logspace(-2, 0, n_depth)
t_rel_list = [0.5, 0.75, 1, 2]
p_rel_list = [1, 10**1, 10**2]

grad_nd_nd, dz_nd = np.meshgrid(grad_nd_nd_list, dz_nd_list, indexing="ij")

grad_nd_gas = grad_nd_nd + 1
spc_work_gas = (grad_nd_gas - 1) * dz_nd
spc_ex_gas = spc_work_gas - np.log((1 + grad_nd_gas * dz_nd)/(1 + dz_nd))
carnot_factor_gas = 1 - 1 / (1 + grad_nd_gas * dz_nd)
ex_eta_gas = spc_ex_gas / (spc_work_gas * carnot_factor_gas)

grad_nd_liq = grad_nd_nd
spc_work_liq = grad_nd_liq * dz_nd
spc_ex_liq = spc_work_liq - np.log(1 + spc_work_liq)
carnot_factor_liq = 1 - 1 / (1 + grad_nd_liq * dz_nd)
ex_eta_liq = spc_ex_liq / (spc_work_liq * carnot_factor_liq)

if np.isnan(ex_eta_gas[0, 0]):

    ex_eta_gas[0, :] = np.ones(ex_eta_gas[0, :].shape)

base_shape = np.array(grad_nd_nd.shape)
res_shape = np.append([len(p_rel_list), len(t_rel_list)], base_shape)

grad_nd_rocks = np.empty(res_shape)
grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
alpha_in = np.empty(res_shape)
t_rocks = np.empty(res_shape)
w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

grad_nd_rocks[:] = np.nan
grad_rocks[:] = np.nan
depth[:] = np.nan
alpha_in[:] = np.nan
t_rocks[:] = np.nan
w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
in_state = PlantThermoPoint(["Methane"], [1], unit_system="MASS BASE SI")
bhe_in = PlantThermoPoint(["Methane"], [1])

states = list()
for i in range(4):
    states.append(in_state.duplicate())

bhe_list = list()
pbar = tqdm(desc="Calculating Points", total=len(p_rel_list) * len(t_rel_list) * len(grad_nd_nd_list) * len(dz_nd_list))
for a in range(len(p_rel_list)):

    p_rel = p_rel_list[a]
    bhe_semi_list = list()

    for i in range(len(t_rel_list)):

        t_in = t_rel_list[i] * in_state.RPHandler.TC
        p_in = p_rel * in_state.RPHandler.PC
        in_state.set_variable("T", t_in)
        in_state.set_variable("P", p_in)
        in_state.copy_state_to(bhe_in)

        dpdt_in = in_state.get_derivative("P", "T", "rho")
        gamma_in = in_state.get_variable("CP/CV")
        v_in = 1 / in_state.get_variable("rho")
        cp_in = in_state.get_variable("cp")
        depth[a, i, :, :] = dz_nd * in_state.get_variable("cp") * t_in / scipy.constants.g

        bhe_sub_list = list()
        for k in range(len(dz_nd_list)):

            bhe = SimplifiedBHE(bhe_in, dz_well=depth[a, i, 0, k], t_rocks=50)
            bhe.update()

            t_input = bhe.points[0].get_variable("T")
            t_down = bhe.points[1].get_variable("T")
            alpha_in[a, i, :, k] = cp_in * (t_down - t_input) / (scipy.constants.g * depth[a, i, 0, k])

            grad_nd_rocks[a, i, :, k] = grad_nd_nd[:, k] + alpha_in[a, i, :, k]
            grad_rocks[a, i, :, k] = grad_nd_rocks[a, i, :, k] * scipy.constants.g / in_state.get_variable("cp")
            t_rocks[a, i, :, k] = t_in + depth[a, i, :, k] * grad_rocks[a, i, :, k]

            bhe_subsub_list = list()
            for j in range(len(grad_nd_nd_list)):

                bhe = SimplifiedBHE(bhe_in, dz_well=depth[a, i, j, k], t_rocks=t_rocks[a, i, j, k]-273.15)
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

                if (not (w_dot < 0 or ex_dot < 0)) and (not states[-1].get_variable("H") < -1e6):

                    cf = 1 - states[0].get_variable("T") / states[2].get_variable("T")

                    w_dot_nds[a, i, j, k] = w_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                    ex_dot_nds[a, i, j, k] = ex_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                    eta_exs[a, i, j, k] = ex_dot / (w_dot * cf)

                    w_dex_mins[a, i, j, k] = (surface_states[1].get_variable("H") - states[0].get_variable("H")) / w_dot
                    w_dex_maxs[a, i, j, k] = (states[3].get_variable("H") - surface_states[1].get_variable("H")) / w_dot

                pbar.update(1)
            bhe_sub_list.append(bhe_subsub_list)
            gc.collect()

        bhe_semi_list.append(bhe_sub_list)
    bhe_list.append(bhe_semi_list)
pbar.close()


# %%-------------------------------------   IDENTIFY INITIAL RISE               -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)
dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
line_styles = ["-", "--", "-.", "-"]
w_dot_rel = w_dot_nds / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

for i in range(len(t_rel_list[:])):

    lines = list()
    for k in range(len(dz_nd_list)):

        lines.append(

            ax.plot(

                grad_nd_rocks[a, i, :, k], w_dot_nds[a, i, :, k], line_styles[i],
                label=dz_label.format(np.log10(dz_nd_list[k])),
                color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

            )[0]

        )

        ax.plot(

            grad_nd_nd[:, k], spc_work_liq[:, k], line_styles[i],
            label=dz_label.format(np.log10(dz_nd_list[k])),
            color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

        )

        ax.plot(

            grad_nd_nd[:, k] + 1, spc_work_gas[:, k], line_styles[i],
            label=dz_label.format(np.log10(dz_nd_list[k])),
            color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

        )

# ax.set_ylim((-np.nanmax(dw_dot_rel)*1.25, np.nanmax(dw_dot_rel)*1.25))
# ax.set_xlim((0.75 * 10 ** -2, 1.25 * 10 ** 4))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("d(${\\dot{w}}^{\\#}$ / ${{\\dot{w}}^{\\#}}_{liq}$) [-]")
ax.set_xlabel("${\\nabla T_{rocks}}^{\\# \\#}$ [-]")
ax.legend(handles=lines, fontsize="12")
plt.show()


# %%-------------------------------------   IDENTIFY DIFFERENT VALUES           -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(13, 5)
cmap = plt.get_cmap('viridis')
t_rel_label = "$T_{{rel}} = 10^{{ {:0.1f} }}$"

w_dot_rel = w_dot_nds / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

lines = list()
min_grad_dz_all = list()
max_grad_dz_all = list()
max_grad_dz_ovr_all = list()

sub_t_rel_list = t_rel_list[:-1]
for i in range(len(sub_t_rel_list)):

    min_grad = list()
    max_grad = list()
    max_grad_ovr = list()

    min_grad_dz = list()
    max_grad_dz = list()
    max_grad_dz_ovr = list()

    for k in range(len(dz_nd_list)):

        max_dw_dot_rel = np.nanmax(dw_dot_rel[a, i, :, k])
        max_w_dot_rel = np.nanmax(w_dot_rel[a, i, np.where(dw_dot_rel[a, i, :, k] > 0), k])
        min_w_dot_rel = np.nanmin(w_dot_rel[a, i, np.where(dw_dot_rel[a, i, :, k] > 0.1 * max_dw_dot_rel), k])

        min_j = np.where(w_dot_rel[a, i, :, k] == min_w_dot_rel)[0][0]
        max_j = np.where(dw_dot_rel[a, i, :, k] == max_dw_dot_rel)[0][0]
        max_ovr_j = np.where(w_dot_rel[a, i, :, k] == max_w_dot_rel)[0][0]

        min_grad.append(grad_nd_nd[min_j, k])
        max_grad.append(grad_nd_nd[max_j, k])
        max_grad_ovr.append(grad_nd_nd[max_ovr_j, k])

        min_grad_dz.append(grad_nd_nd[min_j, k] * dz_nd[min_j, k])
        max_grad_dz.append(grad_nd_nd[max_j, k] * dz_nd[max_j, k])
        max_grad_dz_ovr.append(grad_nd_nd[max_ovr_j, k] * dz_nd[max_ovr_j, k])

    min_grad_dz_all.append(min_grad_dz)
    max_grad_dz_all.append(max_grad_dz)
    max_grad_dz_ovr_all.append(max_grad_dz_ovr)

    lines.append(

        axs[0].plot(

            dz_nd_list, min_grad, line_styles[i],
            label=t_rel_label.format(np.log10(sub_t_rel_list[i])),
            color=cmap(norm(0))

        )[0]

    )

    axs[0].plot(

        dz_nd_list, max_grad, line_styles[i],
        label=t_rel_label.format(np.log10(sub_t_rel_list[i])),
        color=cmap(norm(0.4))

    )

    axs[0].plot(

        dz_nd_list, max_grad_ovr, line_styles[i],
        label=t_rel_label.format(np.log10(sub_t_rel_list[i])),
        color=cmap(norm(0.8))

    )

axs[0].set_xlabel("${{\\Delta z}}^{{\\#}}$ [-]")
axs[0].set_ylabel("${\\nabla T_{rocks}}^{\\# \\#}$ [-]")
axs[0].legend(handles=lines, fontsize="12")
axs[0].set_xscale("log")
axs[0].set_yscale("log")

min_grad_dz_all = np.array(min_grad_dz_all)
max_grad_dz_all = np.array(max_grad_dz_all)
max_grad_dz_ovr_all = np.array(max_grad_dz_ovr_all)

k = 0
axs[1].plot(sub_t_rel_list, min_grad_dz_all[:, k], label="Initial", color=cmap(norm(0)))
axs[1].plot(sub_t_rel_list, max_grad_dz_all[:, k], label="Max Increase", color=cmap(norm(0.4)))
axs[1].plot(sub_t_rel_list, max_grad_dz_ovr_all[:, k], label="Max", color=cmap(norm(0.8)))
axs[1].set_ylabel("${{\\Delta z}}^{{\\#}}{\\nabla T_{rocks}}^{\\#}$ [-]")
axs[1].set_xlabel("$T_{{rel}}$ [-]")
axs[1].set_xscale("log")
axs[1].set_yscale("log")

axs[1].legend(fontsize="12")
plt.tight_layout()
plt.show()
