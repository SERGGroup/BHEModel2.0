# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   DEFINE WORKING FOLDERS                ------------------------------------->
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "1 - Ideal Gas or Liquid System",
    "res"

)


# %%-------------------------------------   RESULTS EVALUATION                  -------------------------------------> #
grad_nd_nd_list = np.logspace(-3, 2, 1000)[1:]
dz_nd_ext = np.logspace(-2, 2, 1000)

dz_nd_list = [10 ** -2, 10 ** -1.5, 10 ** -1, 10 ** -0.5, 1]
grad_nd_nd, dz_nd = np.meshgrid(grad_nd_nd_list, dz_nd_list)

grad_nd_liq = grad_nd_nd
spc_work_liq = grad_nd_liq * dz_nd
spc_ex_liq = spc_work_liq - np.log(1 + spc_work_liq)
carnot_factor = 1 - 1 / (1 + grad_nd_liq * dz_nd)
ex_eta_liq = spc_ex_liq / (spc_work_liq * carnot_factor)

grad_nd_vap = 1 + grad_nd_nd
spc_work_gas = (grad_nd_vap - 1) * dz_nd
spc_ex_gas = spc_work_gas - np.log((1 + grad_nd_vap * dz_nd)/(1 + dz_nd))
carnot_factor = 1 - 1 / (1 + grad_nd_vap * dz_nd)
ex_eta_gas = spc_ex_gas / (spc_work_gas * carnot_factor)

dir_exp_perc_abs_min = dz_nd / (1 + dz_nd*grad_nd_vap)
dir_exp_perc_abs_max = dz_nd/((1 + dz_nd)*(1 + dz_nd*grad_nd_vap - dz_nd))

dir_exp_perc_base = dz_nd_ext / (1 + dz_nd_ext)

dir_exp_perc_rel_min = (1 + dz_nd)/(1 + dz_nd*grad_nd_vap)
dir_exp_perc_rel_max = (1 + dz_nd)/((1 + dz_nd)*(1 + dz_nd*grad_nd_vap - dz_nd))


# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
print_grad_nd_nd = False
colors = [

    "tab:blue", "tab:orange",
    "tab:green", "tab:red",
    "tab:purple"

]

label_str = "${{\\Delta z}}^{{\\#}}\\ =\\ 10^{{ {order} }}$"
fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
# remove the underlying axes
for ax in base_axs[:, 1]:
    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]

lines = [list(), list(), list()]

if print_grad_nd_nd:

    x_values_liq = grad_nd_nd
    x_values_vap = grad_nd_nd

else:

    x_values_liq = grad_nd_liq
    x_values_vap = grad_nd_vap

for i in range(len(dz_nd_list)):

    color = colors[i]
    dz_nd_curr = dz_nd_list[i]
    label_curr = label_str.format(

        order=np.round(np.log10(dz_nd_curr), 1)

    )

    axs[0].plot(x_values_liq[i, :], spc_work_liq[i, :], "--", color=color)
    axs[1].plot(x_values_liq[i, :], spc_ex_liq[i, :], "--", color=color)
    axs[2].plot(x_values_liq[i, :], ex_eta_liq[i, :], "--", color=color)

    lines[0].append(

        axs[0].plot(

            x_values_vap[i, :],
            spc_work_gas[i, :],
            "-", label=label_curr,
            color=color

        )[0]

    )

    lines[1].append(

        axs[1].plot(

            x_values_vap[i, :],
            spc_ex_gas[i, :],
            "-", label=label_curr,
            color=color

        )[0]

    )

    lines[2].append(

        axs[2].plot(

            x_values_vap[i, :],
            ex_eta_gas[i, :],
            "-", label=label_curr,
            color=color

        )[0]

    )

y_names = ["${\\dot{w}}^{\\#}$ [-]", "${\\dot{e}_x}^{\\#}$ [-]", "${\\eta}_{ex}$ [-]"]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])

    if print_grad_nd_nd:

        axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#\\#}$ [-]")

        axs[k].set_xlim((

            np.min(x_values_vap),
            np.max(x_values_vap)

        ))

    else:

        axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#}$ [-]")
        axs[k].set_xlim((

            np.min(x_values_vap)*0.75,
            np.max(x_values_vap)*1.25

        ))

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")

        if print_grad_nd_nd:
            lim_value = spc_ex_liq

        else:
            lim_value = spc_ex_gas

        axs[k].set_ylim((

            np.min(lim_value) / 2,
            np.max(lim_value) * 2

        ))

    if not y_names[k] == y_names[1]:
        axs[k].legend(handles=lines[k], fontsize="8")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   PLOT EXPANSION WORK PERC            -------------------------------------> #
colors = [

    "tab:blue", "tab:orange",
    "tab:green", "tab:red",
    "tab:purple"

]

label_str = "${{\\Delta z}}^{{\\#}}\\ =\\ 10^{{ {order} }}$"
fig, axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 7)
gs = axs[0, 0].get_gridspec()
# remove the underlying axes
for ax in axs[:, 0]:
    ax.remove()

axbig = fig.add_subplot(gs[:, 0])

axs[0, 1].plot(

    dz_nd_ext,
    dir_exp_perc_base,
    "-", color="tab:blue"

)

lines = list()
big_lines = list()
for i in range(len(dz_nd_list)):

    color = colors[i]
    dz_nd_curr = dz_nd_list[i]
    label_curr = label_str.format(

        order=np.round(np.log10(dz_nd_curr), 1)

    )

    axs[1, 1].fill_between(

        x_values_vap[i, :],
        dir_exp_perc_rel_min[i, :],
        dir_exp_perc_rel_max[i, :],
        color=color, alpha=0.1
    )
    axs[1, 1].plot(x_values_vap[i, :], dir_exp_perc_rel_min[i, :], "-", color=color)
    lines.append(

        axs[1, 1].plot(

            x_values_vap[i, :],
            dir_exp_perc_rel_max[i, :],
            "-", label=label_curr,
            color=color

        )[0]

    )

    axbig.fill_between(

        x_values_vap[i, :],
        dir_exp_perc_abs_min[i, :],
        dir_exp_perc_abs_max[i, :],
        color=color, alpha=0.1
    )
    axbig.plot(x_values_vap[i, :], dir_exp_perc_abs_min[i, :], "-", color=color)
    big_lines.append(

        axbig.plot(

            x_values_vap[i, :],
            dir_exp_perc_abs_max[i, :],
            "-", label=label_curr,
            color=color

        )[0]

    )

all_axs = [axbig, axs[0, 1], axs[1, 1]]
y_names = ["$\\dot{w}_{dex}$ [-]", "$\\dot{w}_{{dex}_{base}}$ [-]", "$\\dot{w}_{{dex}_{rel}}$ [-]"]
x_names = ["${\\nabla T_{rocks}}^{\\#} - 1$ [-]", "${\\Delta z}^{\\#}$ [-]", "${\\nabla T_{rocks}}^{\\#} - 1$ [-]"]
x_values = [x_values_vap.ravel(), dz_nd_ext, x_values_vap.ravel()]
for k in range(len(all_axs)):

    all_axs[k].set_xscale("log")
    all_axs[k].set_xlabel(x_names[k])
    all_axs[k].set_ylabel(y_names[k])
    all_axs[k].set_xlim((

        np.min(x_values[k]),
        np.max(x_values[k])

    ))

axbig.legend(handles=big_lines, fontsize="10")
axs[1, 1].legend(handles=lines, fontsize="10")

props = dict(boxstyle='round', facecolor='white', alpha=0.75)
axbig.text(

    0.9, 0.65,
    '$\\dot{w}_{dex} = \\dot{w}_{{dex}_{base}} \\dot{w}_{{dex}_{rel}}$',
    verticalalignment='top', horizontalalignment="right", bbox=props,
    fontsize = 12, transform=axbig.transAxes
)
plt.tight_layout(pad=2)
plt.show()
