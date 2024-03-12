# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.signal import argrelextrema
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import pandas as pd
import numpy as np


# %%-------------------------------------   DEFINE CONDITIONS                   -------------------------------------> #
fluid = "Ethane"
n_grad = 350
n_depth = 20
n_t_rel = 3

grad_nd_nd_list = np.logspace(-2, 4.5, n_grad)
dz_nd_list = np.logspace(-5, -1, n_depth)
t_rel_list = [0.4, 0.5, 0.75]


# %%-------------------------------------   EVALUATE COEFFICIENTS               -------------------------------------> #
p_sat_list = np.logspace(-3, 0, 40000)
pbar = tqdm(desc="Calculating Saturation Points", total=len(p_sat_list))
t_sat_arr = np.empty(np.array([1, len(p_sat_list)]))
t_sat_arr[:, :] = np.nan

tp_vap = PlantThermoPoint([fluid], [1], unit_system="MASS BASE SI")
tp_liq = tp_vap.duplicate()

p_crit = tp_vap.RPHandler.PC
t_crit = tp_vap.RPHandler.TC

for j in range(len(p_sat_list)):

    tp_vap.set_variable("P", p_sat_list[j]*p_crit)
    tp_vap.set_variable("Q", 1)
    tp_liq.set_variable("P", p_sat_list[j] * p_crit)
    tp_liq.set_variable("Q", 0)

    if tp_vap.get_variable("T") / t_crit > 0.4:
        t_sat_arr[0, j] = tp_vap.get_variable("T") / t_crit

    pbar.update(1)

pbar.close()
acntr_factors = [tp_vap.evaluate_RP_code("ACF")]

def opt_function(i_curr, minimise_max=False):

    def funct(x):

        a_curr = x[0]
        b_curr = np.log(0.7) * a_curr / (np.power(10, - a_curr * (acntr_factors[i_curr] + 1)) - 1)
        t_sat_actr = np.exp(b_curr / a_curr * (np.power(p_sat_list, a_curr) - 1))

        if minimise_max:
            return np.nanmax(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

        return np.nanmean(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

    return funct

res = minimize(opt_function(0), np.array([0.11]))
a_fluid = res.x[0]
b_fluid = np.log(0.7) * a_fluid / (np.power(10, - a_fluid * (acntr_factors[0] + 1)) - 1)


# %%-------------------------------------   INIT POINTS                         -------------------------------------> #
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
betas_des = np.empty(res_shape)
w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

grad_nd_rocks[:] = np.nan
grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan
betas_des[:] = np.nan
w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE CORRELATION COEFFICIENTS  -------------------------------------> #
in_state = PlantThermoPoint([fluid], [1], unit_system="MASS BASE SI")
calc_coefficients = list()
p_in_list = list()
dz_coeff = list()

pbar = tqdm(desc="Calculating Correlation Coefficients", total=len(t_rel_list))

for i in range(len(t_rel_list)):

    t_in = t_rel_list[i] * in_state.RPHandler.TC
    in_state.set_variable("T", t_in)
    in_state.set_variable("Q", 0)
    p_in = in_state.get_variable("P") * 1.0001

    in_state.set_variable("T", t_in)
    in_state.set_variable("P", p_in)

    dpdt_in = in_state.get_derivative("P", "T", "rho")
    gamma_in = in_state.get_variable("CP/CV")
    v_in = 1 / in_state.get_variable("rho")
    cp_in = in_state.get_variable("cp")

    rho_in = in_state.get_variable("rho")
    dpdrho_t = in_state.get_derivative("P", "rho", "T")

    dpdv_t = - rho_in ** 2 * dpdrho_t
    dvdp = 1 / (gamma_in * dpdv_t)

    alpha_nd = cp_in * rho_in * t_in / p_in
    beta_nd = (p_in * rho_in * dvdp + 1) / a_fluid

    p_in_list.append(p_in)
    dz_coeff.append(rho_in * cp_in * t_in / p_in)
    calc_coefficients.append([

        alpha_nd * b_fluid * np.power(p_in / in_state.RPHandler.PC, a_fluid),
        a_fluid * alpha_nd * (1 - beta_nd),
        3 * a_fluid * alpha_nd * beta_nd / (beta_nd - 1),
        - a_fluid * alpha_nd * beta_nd

    ])

    pbar.update(1)

pbar.close()


# %%-------------------------------------   CALCULATE POINTS                    -------------------------------------> #
bhe_in = PlantThermoPoint([fluid], [1])

states = list()

for i in range(4):
    states.append(in_state.duplicate())

bhe_list = list()
pbar = tqdm(desc="Calculating Points", total=len(t_rel_list) * len(grad_nd_nd_list) * len(dz_nd_list))
for i in range(len(t_rel_list)):

    t_in = t_rel_list[i] * in_state.RPHandler.TC
    in_state.set_variable("T", t_in)
    in_state.set_variable("Q", 0)
    p_in = in_state.get_variable("P") * 1.0001

    in_state.set_variable("T", t_in)
    in_state.set_variable("P", p_in)
    in_state.copy_state_to(bhe_in)

    dpdt_in = in_state.get_derivative("P", "T", "rho")
    gamma_in = in_state.get_variable("CP/CV")
    v_in = 1 / in_state.get_variable("rho")
    cp_in = in_state.get_variable("cp")

    depth[i, :, :] = dz_nd * in_state.get_variable("cp") * t_in / scipy.constants.g

    bhe_sub_list = list()
    for k in range(len(dz_nd_list)):

        bhe = SimplifiedBHE(bhe_in, dz_well=depth[i, 0, k], t_rocks=50)
        bhe.update()

        alpha_in = cp_in * (bhe.points[1].get_variable("T") - bhe.points[0].get_variable("T")) / (scipy.constants.g * depth[i, 0, k])
        grad_nd_rocks[i, :, k] = grad_nd_nd[:, k] + alpha_in
        grad_rocks[i, :, k] = grad_nd_rocks[i, :, k] * scipy.constants.g / in_state.get_variable("cp")
        t_rocks[i, :, k] = t_in + depth[i, :, k] * grad_rocks[i, :, k]

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
            cf = 1 - states[0].get_variable("T") / states[2].get_variable("T")

            betas_des[i, j, k] = bhe.points[1].get_variable("P") / bhe.points[0].get_variable("P")
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

        wrong_indices = np.where(eta_exs[i, :, k] > 1)
        w_dot_nds[i, wrong_indices, k] = np.nan
        ex_dot_nds[i, wrong_indices, k] = np.nan
        eta_exs[i, wrong_indices, k] = np.nan
        w_dex_mins[i, wrong_indices, k] = np.nan
        w_dex_maxs[i, wrong_indices, k] = np.nan



# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
i = 0

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

                x_values = grad_nd_rocks[i, :, :]

                lines[m].append(

                    axs[m].plot(

                        x_values[:, k], cmp_y_values[n][m][i, :, k], "-",
                        label=dz_label.format(np.log10(dz_nd_list[k])),
                        color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

                    )[0]

                )

        else:

            if cmp_labels[n] == "Ideal Gas":
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
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#}$ [-]")

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")

    else:
        axs[k].legend(handles=lines[k], fontsize="8")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   PLOT DIFFERENCE WITH LIQUID V2.0    -------------------------------------> #
i = 0
m = 0

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
x_values = grad_nd_rocks[i, :, :]
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

        x_values[:, k] * dz_nd[:, k] + 1, cmp_y_values[0][m][i, :, k] / liq_value[:, k], "-",
        label=dz_label.format(np.log10(dz_nd_list[k])),
        color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

    )


y_names = [

    "${\\dot{w}}^{\\#}$ / ${{\\dot{w}}^{\\#}}_{liq}$ [-]",
    "${\\dot{e}_x}^{\\#}$ / ${{\\dot{e}_x}^{\\#}}_{liq}$ [-]",
    "${\\eta}_{ex}$ / ${{\\eta}_{ex}}_{liq}$ [-]"

]

axs[0].legend(handles=lines[m], fontsize="8")
base_x_label = "${\\nabla T_{rocks}}^{\\#}$"

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[m])
    axs[k].set_xlabel("{} [-]".format(base_x_label))
    base_x_label += "${\\Delta z}^{\\#}$"


plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   IDENTIFY INITIAL RISE V2.0          -------------------------------------> #
i = 0
m = 0
mean_window = 1

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
x_values = grad_nd_rocks[i, :, :]
grad_nd_liq = x_values

liq_value = grad_nd_liq * dz_nd
j_min = np.min(np.where(grad_nd_nd[:, 0] > 1))
liq_difference = cmp_y_values[0][m][i, :, :] / liq_value[:, :]

for k in range(len(dz_nd_list)):

    liq_difference_mean = pd.Series(liq_difference[:, k]).rolling(mean_window, center=True).mean().iloc[:].values
    delta_liq_diff = liq_difference_mean[1:] - liq_difference_mean[:-1]
    delta_liq_diff_var = (k/2+1)*(np.nanmax(delta_liq_diff) - np.nanmin(delta_liq_diff))
    delta_liq_diff_norm = delta_liq_diff / delta_liq_diff_var
    print(argrelextrema(delta_liq_diff_norm, np.greater))
    lines[m].append(

        axs[0].plot(

            x_values[:-1, k], delta_liq_diff_norm, "-",
            label=dz_label.format(np.log10(dz_nd_list[k])),
            color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

        )[0]

    )

    axs[1].plot(

        x_values[:-1, k] * dz_nd[:-1, k] + 1, delta_liq_diff_norm, "-",
        label=dz_label.format(np.log10(dz_nd_list[k])),
        color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

    )

min_values = calc_coefficients[i][0] * dz_nd[0, :]
y_values = np.zeros(min_values.shape)

y_names = [

    "{\\dot{w}}^{\\#} / {{\\dot{w}}^{\\#}}_{liq}",
    "{\\dot{e}_x}^{\\#} / {{\\dot{e}_x}^{\\#}}_{liq}",
    "{\\eta}_{ex} / {{\\eta}_{ex}}_{liq}"

]

derivative_append = "$d({y_value}) / d({x_value})$"

axs[0].legend(handles=lines[m], fontsize="8")
base_x_label = "{\\nabla T_{rocks}}^{\\#}"

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(derivative_append.format(y_value=y_names[m], x_value=base_x_label))
    axs[k].set_xlabel("{} [-]".format("${}$".format(base_x_label)))
    base_x_label += "{\\Delta z}^{\\#}"


plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   IDENTIFY DIFFERENT VALUES           -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(7, 5)
cmap = plt.get_cmap('viridis')
t_rel_label = "$T_{{rel}} = 10^{{ {:0.1f} }}$"

w_dot_rel = w_dot_nds / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

lines = list()
max_grad_dz_all = list()
line_styles = ["-", "--", "-."]

sub_t_rel_list = t_rel_list
for i in range(len(sub_t_rel_list)):

    max_grad = list()

    for k in range(len(dz_nd_list) - 2):

        max_dw_dot_rel = np.nanmax(dw_dot_rel[i, :, k])
        max_j = np.where(dw_dot_rel[i, :, k] == max_dw_dot_rel)[0][0]
        max_grad.append(grad_nd_rocks[i, max_j, k] * dz_nd[max_j, k] + 1)

    max_grad_dz_all.append(max_grad)
    max_grad = np.array(max_grad)

    coeff = calc_coefficients[i]
    approx = coeff[0] * dz_nd_list[: -2]

    for n in range(0):
        approx = approx * (1 + coeff[n + 1] * dz_nd_list[: -2] / (n + 2))

    # line = ax.plot(np.log(dz_nd_list[:-2] + 1), np.log(max_grad), label=t_rel_label.format(sub_t_rel_list[i]))[0]
    # ax.plot(np.log(dz_nd_list[:-2] + 1), approx, "--", color=line.get_color())

    # x_values = np.power(betas_des[i, 0, :-2], a_fluid) - 1
    x_values = np.power(np.log(dz_nd_list[:-2] + 1) * dz_coeff[i]*(dz_nd_list[:-2] + 1) + 1, a_fluid) - 1
    y_values = np.log(max_grad) * a_fluid / (b_fluid * np.power(p_in_list[i], a_fluid))
    ax.plot(x_values, y_values, label=t_rel_label.format(sub_t_rel_list[i]))

    # ax.plot(np.log(dz_nd_list+1) * (dz_coeff[i]*(dz_nd_list + 1)) + 1, betas_des[i, 0, :], label=t_rel_label.format(sub_t_rel_list[i]))

# ax.set_ylim([0, 10])
# ax.set_xlim([0, 10])
plt.legend()
plt.tight_layout()
plt.show()


