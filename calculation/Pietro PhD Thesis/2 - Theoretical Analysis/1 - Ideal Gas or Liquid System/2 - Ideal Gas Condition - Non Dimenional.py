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


# %%-------------------------------------   IDEAL GAS EVALUATION                -------------------------------------> #
n_depths = 5
n_grads = 1000
depths_nd_list = np.logspace(-2, 0, n_depths)
grads_nd_list = np.logspace(0, 3, n_grads + 1)[1:]
grads_nd, depths_nd = np.meshgrid(grads_nd_list, depths_nd_list, indexing='ij')
gamma = 1.4             # -

grads_nd_new = grads_nd * depths_nd

t_2 = 1 + grads_nd_new
tau_des = 1 + depths_nd
beta_des = tau_des ** (gamma / (gamma - 1))

tau_asc = (1 - depths_nd / t_2)
beta_asc = tau_asc ** (gamma / (gamma - 1))
t_3 = t_2 * tau_asc

beta_turb = beta_asc * beta_des
tau_turb = beta_turb ** ((gamma - 1) / gamma)

t_3ii = t_3 / tau_turb
t_3i = tau_turb

w_exp_nd_max_exp = (t_3 - t_3ii)
w_exp_nd_min_exp = (t_3i - 1)
q_cool_nd_max_exp = (t_3 - t_3i)
q_cool_nd_min_exp = (t_3ii - 1)

# Which can be calculated directly considering that:
t_3ii = (1 + grads_nd * depths_nd) / (1 + depths_nd)
t_3i = (1 + grads_nd * depths_nd - depths_nd) * (1 + depths_nd) / (1 + grads_nd * depths_nd)

# Resulting in:
grads_nd_rel = grads_nd - 1
w_exp_nd_max = grads_nd_rel * depths_nd ** 2 / (1 + depths_nd)
w_exp_nd_min = grads_nd_rel * depths_nd ** 2 / (1 + grads_nd * depths_nd)
q_cool_nd_max = grads_nd_rel * depths_nd / (1 + grads_nd * depths_nd) * (1 + grads_nd * depths_nd - depths_nd)
q_cool_nd_min = grads_nd_rel * depths_nd / (1 + depths_nd)


# %%-------------------------------------   PLOT NON-DIMENSIONAL RESULT         -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=400)
fig.set_size_inches(10, 4)
label_format = "${{\\Delta z}}_{{nd}} = 10^{{{:0.1f}}}$"
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
width = 1.5
alpha = 0.2

fill_line = [list(), list()]
for i in range(len(depths_nd_list)):

    axs[0].fill_between(grads_nd_rel[:, i], w_exp_nd_max[:, i], w_exp_nd_min[:, i], color=colors[i], alpha=alpha)
    axs[0].plot(grads_nd_rel[:, i], w_exp_nd_max[:, i], color=colors[i], linewidth=width)
    fill_line[0].append(

        axs[0].plot(

            grads_nd_rel[:, i], w_exp_nd_min[:, i],
            label=label_format.format(np.log10(depths_nd_list[i])),
            color=colors[i], linewidth=width

        )[0]

    )

    axs[1].fill_between(grads_nd_rel[:, i], q_cool_nd_max[:, i], q_cool_nd_min[:, i], color=colors[i], alpha=alpha)
    axs[1].plot(grads_nd_rel[:, i], q_cool_nd_max[:, i], color=colors[i], linewidth=width)
    fill_line[1].append(

        axs[1].plot(

            grads_nd_rel[:, i], q_cool_nd_min[:, i],
            label=label_format.format(np.log10(depths_nd_list[i])),
            color=colors[i], linewidth=width

        )[0]

    )

for j in range(2):

    axs[j].set_xscale("log")
    axs[j].set_yscale("log")

    axs[j].set_xlabel("${\\nabla T_{rocks}}_{\\ nd}$ - 1 [-]")
    axs[j].set_xlim(np.min(grads_nd_rel), np.max(grads_nd_rel))

    labs = [l.get_label() for l in fill_line[j]]
    axs[j].legend(fill_line[j], labs, fontsize=8)


axs[0].set_title("Expansion Power")
axs[0].set_ylabel("${\\dot{w}}_{nd}$ [-]")

axs[1].set_title("Cooling Power")
axs[1].set_ylabel("${\\dot{q}}_{nd}$ [-]")

plt.tight_layout()
plt.show()


# %%-------------------------------------   EVALUATE MODIFIERS                  -------------------------------------> #
w_exp_base = grads_nd_rel * depths_nd ** 2 / (1 + depths_nd)
q_cool_base = grads_nd_rel * depths_nd / (1 + depths_nd)

w_exp_mod_max = np.ones(depths_nd.shape)
w_exp_mod_min = (1 + depths_nd) / (1 + grads_nd * depths_nd)
q_cool_mod_max = ((1 + depths_nd) * (1 + grads_nd * depths_nd - depths_nd) / (1 + grads_nd * depths_nd) - 1) / depths_nd
q_cool_mod_min = np.zeros(depths_nd.shape)


# %%-------------------------------------   PLOT MODIFIERS                      -------------------------------------> #
fig, axs = plt.subplots(2, 2, dpi=400)
fig.set_size_inches(10, 6)
label_format = "${{\\Delta z}}_{{nd}} = 10^{{{:0.1f}}}$"
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
width = 1.5
alpha = 0.1

fill_line = [list(), list()]

for i in range(len(depths_nd_list)):

    axs[0][0].plot(grads_nd_rel[:, i] , w_exp_base[:, i], color=colors[i], linewidth=width)
    axs[0][1].plot(grads_nd_rel[:, i], q_cool_base[:, i], color=colors[i], linewidth=width)

    depth_curr = depths_nd_list[i]
    log_depth = np.log10(depth_curr)

    axs[1][0].fill_between(grads_nd_rel[:, i], w_exp_mod_max[:, i], w_exp_mod_min[:, i], color=colors[i], alpha=alpha)
    axs[1][0].plot(grads_nd_rel[:, i], w_exp_mod_max[:, i], color=colors[i], linewidth=width)
    fill_line[0].append(

        axs[1][0].plot(

            grads_nd_rel[:, i], w_exp_mod_min[:, i],
            label=label_format.format(log_depth),
            color=colors[i], linewidth=width

        )[0]

    )

    axs[1][1].fill_between(grads_nd_rel[:, i], q_cool_mod_max[:, i], q_cool_mod_min[:, i], color=colors[i], alpha=alpha)
    axs[1][1].plot(grads_nd_rel[:, i], q_cool_mod_max[:, i], color=colors[i], linewidth=width)
    fill_line[1].append(

        axs[1][1].plot(

            grads_nd_rel[:, i], q_cool_mod_min[:, i],
            label=label_format.format(log_depth),
            color=colors[i], linewidth=width

        )[0]

    )

label_names = ["w", "q"]
modifier_name = ["base", "mod"]

axs[0][0].set_title("Expansion Power")
axs[0][1].set_title("Cooling Power")

for j in range(2):

    axs[0][j].set_yscale("log")
    # axs[1][j].set_yscale("log")

    for k in range(2):

        axs[k][j].set_xscale("log")
        axs[k][j].set_xlim(np.min(grads_nd_rel), np.max(grads_nd_rel))

        axs[k][j].set_xlabel("${\\nabla T_{rocks}}_{\\ nd}$ - 1 [-]")
        axs[k][j].set_ylabel("${{ \\dot{{ {name} }}_{{ {modifier} }} }}$ [-]".format(

            name=label_names[j],
            modifier=modifier_name[k]

        ))

        labs = [l.get_label() for l in fill_line[j]]
        axs[k][j].legend(fill_line[j], labs, fontsize=8)

plt.tight_layout()
plt.show()

