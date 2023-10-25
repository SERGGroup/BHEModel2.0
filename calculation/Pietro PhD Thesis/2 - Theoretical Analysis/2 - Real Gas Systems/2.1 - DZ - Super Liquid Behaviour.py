# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 350
n_depth = 5
n_t_rel = 3

grad_nd_nd_list = np.logspace(-2, 4.5, n_grad)
dz_nd_list = np.logspace(-2, 0, n_depth)
t_rel_list = [0.5, 0.75, 1, 2]

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
res_shape = np.append([len(t_rel_list)], base_shape)

grad_nd_rocks = np.empty(res_shape)
grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
t_rocks = np.empty(res_shape)
w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

grad_nd_rocks[:] = np.nan
grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan
w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
p_rel = 2*10**1

in_state = PlantThermoPoint(["Ethane"], [1], unit_system="MASS BASE SI")
bhe_in = PlantThermoPoint(["Ethane"], [1])

states = list()
for i in range(4):
    states.append(in_state.duplicate())

bhe_list = list()
pbar = tqdm(desc="Calculating Points", total=len(t_rel_list) * len(grad_nd_nd_list) * len(dz_nd_list))
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
    alpha_in = cp_in * (1 - 1 / gamma_in) / (v_in * dpdt_in)

    depth[i, :, :] = dz_nd * in_state.get_variable("cp") * t_in / scipy.constants.g
    t_rocks[i, :, :] = t_in + depth[i, :, :] * grad_rocks[i, :, :]

    bhe_sub_list = list()
    for k in range(len(dz_nd_list)):

        bhe = SimplifiedBHE(bhe_in, dz_well=depth[i, 0, k], t_rocks=50)
        bhe.update()

        alpha_in = cp_in * (bhe.points[1].get_variable("T") - bhe.points[0].get_variable("T")) / (scipy.constants.g * depth[i, 0, k])
        grad_nd_rocks[i, :, k] = grad_nd_nd[:, k] + alpha_in
        grad_rocks[i, :, k] = grad_nd_rocks[i, :, k] * scipy.constants.g / in_state.get_variable("cp")

        bhe_subsub_list = list()
        for j in range(len(grad_nd_nd_list)):

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

            if True: # (not (w_dot < 0 or ex_dot < 0)) and (not states[-1].get_variable("H") < -1e6):

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

# %%-------------------------------------   REMOVE WRONG VALUES                 -------------------------------------> #
for i in range(len(t_rel_list)):

    for k in range(len(dz_nd_list)):

        wrong_indices = np.where(eta_exs[i, :, k] < 0)
        w_dot_nds[i, wrong_indices, k] = np.nan
        ex_dot_nds[i, wrong_indices, k] = np.nan
        eta_exs[i, wrong_indices, k] = np.nan
        w_dex_mins[i, wrong_indices, k] = np.nan
        w_dex_maxs[i, wrong_indices, k] = np.nan


# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
i = 0
x_grad_nd = True

fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

for ax in base_axs[:, 1]:

    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]

dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
cmp_labels = ["Liquid", dz_label, "Ideal Gas"]
cmp_y_values = [

    [spc_work_liq, spc_ex_liq, ex_eta_liq],
    [w_dot_nds, ex_dot_nds, eta_exs],
    [spc_work_gas, spc_ex_gas, ex_eta_gas],

]

print(t_rel_list[i])
lines = [list(), list(), list()]
for m in range(len(cmp_y_values[0])):

    for n in range(len(cmp_labels)):

        if cmp_labels[n] == dz_label:

            for k in range(len(dz_nd_list)):

                if x_grad_nd:
                    x_values = grad_nd_rocks[i, :, :]

                else:
                    x_values = grad_nd_nd

                lines[m].append(

                    axs[m].plot(

                        x_values[:, k], cmp_y_values[n][m][i, :, k], "-",
                        label=dz_label.format(np.log10(dz_nd_list[k])),
                        color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

                    )[0]

                )

        else:

            if x_grad_nd and cmp_labels[n] == "Ideal Gas":
                x_values = grad_nd_nd + 1

            else:
                x_values = grad_nd_nd

            for k in range(len(dz_nd_list)):

                line = axs[m].plot(

                    x_values[:, k], cmp_y_values[n][m][:, k], "-",
                    label=cmp_labels[n], color=cmap(norm(n)),
                    linewidth=2

                )[0]

                if k == 0:
                    lines[m].append(line)


y_names = ["${\\dot{w}}^{\\#}$ [-]", "${\\dot{e}_x}^{\\#}$ [-]", "${\\eta}_{ex}$ [-]"]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\# \\#}$ [-]")

    if x_grad_nd:
        axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#}$ [-]")

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")

    else:
        axs[k].legend(handles=lines[k], fontsize="8")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   PLOT DIFFERENCE WITH LIQUID         -------------------------------------> #
i = 0
x_grad_nd = False

fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

for ax in base_axs[:, 1]:

    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]

dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
cmp_y_values = [

    [w_dot_nds, ex_dot_nds, eta_exs],
    [spc_work_liq, spc_ex_liq, ex_eta_liq],
    [spc_work_gas, spc_ex_gas, ex_eta_gas]

]
print(t_rel_list[i])
lines = [list(), list(), list()]
for m in range(len(cmp_y_values[0])):

    if x_grad_nd:
        x_values = grad_nd_rocks[i, :, :]
        grad_nd_liq = x_values
    else:
        x_values = grad_nd_nd

    grad_nd_liq = x_values

    if m == 0:
        liq_value = grad_nd_liq * dz_nd

    elif m == 1:
        liq_value = spc_work_liq - np.log(1 + spc_work_liq)

    else:
        carnot_factor_liq = 1 - 1 / (1 + grad_nd_liq * dz_nd)
        liq_value = spc_ex_liq / (spc_work_liq * carnot_factor_liq)



    for k in range(len(dz_nd_list)):

        lines[m].append(

            axs[m].plot(

                x_values[:, k], cmp_y_values[0][m][i, :, k] / liq_value[:, k], "-",
                label=dz_label.format(np.log10(dz_nd_list[k])),
                color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

            )[0]

        )

    if x_grad_nd:
        x_values = grad_nd_nd + 1

    else:
        x_values = grad_nd_nd

    grad_nd_liq = x_values

    if m == 0:
        liq_value = grad_nd_liq * dz_nd

    elif m == 1:
        liq_value = spc_work_liq - np.log(1 + spc_work_liq)

    else:
        carnot_factor_liq = 1 - 1 / (1 + grad_nd_liq * dz_nd)
        liq_value = spc_ex_liq / (spc_work_liq * carnot_factor_liq)

    line = axs[m].plot(

        x_values[:, 0], cmp_y_values[2][m][:, 0] / liq_value[:, 0], "-",
        label="Ideal Gas", color=cmap(norm(1)),
        linewidth=2

    )[0]

    lines[m].append(line)

y_names = [

    "${\\dot{w}}^{\\#}$ / ${{\\dot{w}}^{\\#}}_{liq}$ [-]",
    "${\\dot{e}_x}^{\\#}$ / ${{\\dot{e}_x}^{\\#}}_{liq}$ [-]",
    "${\\eta}_{ex}$ / ${{\\eta}_{ex}}_{liq}$ [-]"

]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\# \\#}$ [-]")

    if x_grad_nd:
        axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#}$ [-]")

    if y_names[k] == y_names[-1]:
        axs[k].legend(handles=lines[k], fontsize="8")

# axs[0].set_ylim(0.72, 2.35)
axs[1].set_yscale("log")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()

# %%-------------------------------------   PLOT DIFFERENCE WITH LIQUID V2.0    -------------------------------------> #
i = 0
m = 0
x_grad_nd = True

fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 4)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
cmp_y_values = [

    [w_dot_nds, ex_dot_nds, eta_exs],
    [spc_work_liq, spc_ex_liq, ex_eta_liq],
    [spc_work_gas, spc_ex_gas, ex_eta_gas]

]
print(t_rel_list[i])
lines = [list(), list(), list()]
if x_grad_nd:
    x_values = grad_nd_rocks[i, :, :]
    grad_nd_liq = x_values
else:
    x_values = grad_nd_nd

grad_nd_liq = x_values

if m == 0:
    liq_value = grad_nd_liq * dz_nd

elif m == 1:
    liq_value = spc_work_liq - np.log(1 + spc_work_liq)

else:
    carnot_factor_liq = 1 - 1 / (1 + grad_nd_liq * dz_nd)
    liq_value = spc_ex_liq / (spc_work_liq * carnot_factor_liq)

for k in range(len(dz_nd_list)):
    lines[m].append(

        axs[0].plot(

            x_values[:, k], cmp_y_values[0][m][i, :, k] / liq_value[:, k], "-",
            label=dz_label.format(np.log10(dz_nd_list[k])),
            color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

        )[0]

    )

    axs[1].plot(

        x_values[:, k] * dz_nd[:, k], cmp_y_values[0][m][i, :, k] / liq_value[:, k], "-",
        label=dz_label.format(np.log10(dz_nd_list[k])),
        color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

    )


y_names = [

    "${\\dot{w}}^{\\#}$ / ${{\\dot{w}}^{\\#}}_{liq}$ [-]",
    "${\\dot{e}_x}^{\\#}$ / ${{\\dot{e}_x}^{\\#}}_{liq}$ [-]",
    "${\\eta}_{ex}$ / ${{\\eta}_{ex}}_{liq}$ [-]"

]

axs[0].legend(handles=lines[m], fontsize="8")
base_x_label = "${\\nabla T_{rocks}}^{\\# \\#}$"
if x_grad_nd:
    base_x_label = "${\\nabla T_{rocks}}^{\\#}$"

for k in range(len(axs)):

    axs[k].set_xscale("log")
    # axs[k].set_yscale("log")
    axs[k].set_ylabel(y_names[m])
    axs[k].set_xlabel("{} [-]".format(base_x_label))
    base_x_label += "${\\Delta z}^{\\#}$"


plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()

# %%-------------------------------------   IDENTIFY INITIAL RISE               -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)
dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
line_styles = ["-", "--", "-."]
w_dot_rel = w_dot_nds / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

for i in range(len(t_rel_list[:])):

    lines = list()
    for k in range(len(dz_nd_list)):
        lines.append(

            ax.plot(

                grad_nd_nd[:-1, k], dw_dot_rel[i, :, k], line_styles[i],
                label=dz_label.format(np.log10(dz_nd_list[k])),
                color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

            )[0]

        )

ax.set_ylim((-np.nanmax(dw_dot_rel)*1.25, np.nanmax(dw_dot_rel)*1.25))
ax.set_xlim((0.75 * 10 ** -2, 1.25 * 10 ** 4))
ax.set_xscale("log")
ax.set_ylabel("d(${\\dot{w}}^{\\#}$ / ${{\\dot{w}}^{\\#}}_{liq}$) [-]")
ax.set_xlabel("${\\nabla T_{rocks}}^{\\# \\#}$ [-]")
ax.legend(handles=lines, fontsize="12")
ax.set_xscale("log")
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

        max_dw_dot_rel = np.nanmax(dw_dot_rel[i, :, k])
        max_w_dot_rel = np.nanmax(w_dot_rel[i, np.where(dw_dot_rel[i, :, k] > 0), k])
        min_w_dot_rel = np.nanmin(w_dot_rel[i, np.where(dw_dot_rel[i, :, k] > 0.1 * max_dw_dot_rel), k])

        min_j = np.where(w_dot_rel[i, :, k] == min_w_dot_rel)[0][0]
        max_j = np.where(dw_dot_rel[i, :, k] == max_dw_dot_rel)[0][0]
        max_ovr_j = np.where(w_dot_rel[i, :, k] == max_w_dot_rel)[0][0]

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
