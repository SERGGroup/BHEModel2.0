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
g = 9.81                # m / s ^ 2
cp = 1000               # J / (kg K)
gamma = 1.4             # -
t_amb = 20 + 273.15     # K

depth_list = np.logspace(2, 4, 100)
grad_list = np.logspace(1, 2, 5) / 1000
depths, grads = np.meshgrid(depth_list, grad_list, indexing='ij')

t_2 = t_amb + grads * depths

tau_des = 1 + g * depths / (cp * t_amb)
beta_des = tau_des ** (gamma / (gamma - 1))

tau_asc = (1 - g * depths / (cp * t_2))
beta_asc = tau_asc ** (gamma / (gamma - 1))
t_3 = t_2 * tau_asc

beta_turb = beta_asc * beta_des
tau_turb = beta_turb ** ((gamma - 1) / gamma)

t_3ii = t_3 / tau_turb
t_3i = t_amb * tau_turb

w_exp_max = cp * (t_3 - t_3ii) / 1000
w_exp_min = cp * (t_3i - t_amb) / 1000
q_cool_max = cp * (t_3 - t_3i) / 1000
q_cool_min = cp * (t_3ii - t_amb) / 1000


# %%-------------------------------------   PLOT DIMENSIONAL RESULT             -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=400)
fig.set_size_inches(10, 4)
label_format = "$\\nabla T_{{rocks}} = {:0.0f}°C/km$"
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
width = 1.5
alpha = 0.2

fill_line = [list(), list()]
for i in range(len(grad_list)):

    axs[0].fill_between(depths[:, i], w_exp_max[:, i], w_exp_min[:, i], color=colors[i], alpha=alpha)
    axs[0].plot(depths[:, i], w_exp_max[:, i], color=colors[i], linewidth=width)
    fill_line[0].append(

        axs[0].plot(

            depths[:, i], w_exp_min[:, i],
            label=label_format.format(grad_list[i] * 1000),
            color=colors[i], linewidth=width

        )[0]

    )

    axs[1].fill_between(depths[:, i], q_cool_max[:, i], q_cool_min[:, i], color=colors[i], alpha=alpha)
    axs[1].plot(depths[:, i], q_cool_max[:, i], color=colors[i], linewidth=width)
    fill_line[1].append(

        axs[1].plot(

            depths[:, i], q_cool_min[:, i],
            label=label_format.format(grad_list[i] * 1000),
            color=colors[i], linewidth=width

        )[0]

    )

for j in range(2):

    axs[j].set_xscale("log")
    axs[j].set_yscale("log")

    axs[j].set_xlabel("$\\Delta z$ [m]")
    axs[j].set_xlim(np.min(depths), np.max(depths))

    labs = [l.get_label() for l in fill_line[j]]
    axs[j].legend(fill_line[j], labs, fontsize=8)

axs[0].set_title("Expansion Power")
axs[0].set_ylabel("$\\dot{w}$ [kJ/kg]")

axs[1].set_title("Cooling Power")
axs[1].set_ylabel("$\\dot{q}$ [kJ/kg]")

plt.tight_layout()
plt.show()
