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
grad_nd_nd_list = np.concatenate(([0], np.logspace(-3, 1, 1000)))
dz_nd_list = [10 ** -1]
grad_nd_nd, dz_nd = np.meshgrid(grad_nd_nd_list, dz_nd_list)

grad_nd_liq = grad_nd_nd
t_nd_liq = 1 + grad_nd_liq * dz_nd
s_nd_liq = np.log(t_nd_liq)
t_nd_liq = t_nd_liq - 1

grad_nd_gas = 1 + grad_nd_nd
t_nd_gas = 1 + grad_nd_gas * dz_nd
s_nd_gas = np.log(t_nd_gas / (1 + dz_nd))
t_3_gas = grad_nd_gas * dz_nd - dz_nd
t_nd_gas = t_nd_gas-1


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(8, 4)
n_lims = [500, -1]

for i in range(2):

    n_lim = n_lims[i]
    t_nd_liq_rel = t_nd_liq / t_nd_liq[0, n_lim]
    t_nd_gas_rel = t_nd_gas / t_nd_gas[0, n_lim]
    t_3_gas_rel = t_3_gas / t_nd_gas[0, n_lim]

    s_liq_rel = s_nd_liq / s_nd_liq[0, n_lim]
    s_gas_rel = s_nd_gas / s_nd_gas[0, n_lim]

    axs[i].plot(s_gas_rel[0, :n_lim], t_nd_gas_rel[0, :n_lim], color="tab:orange", label="Ideal Gas", zorder=1)
    axs[i].plot(s_liq_rel[0, :n_lim], t_3_gas_rel[0, :n_lim], "--", color="tab:orange", zorder=1)
    axs[i].plot(s_liq_rel[0, :n_lim], t_nd_liq_rel[0, :n_lim], color="tab:blue", label="Liquid", zorder=1)

    axs[i].plot([0, 0], [0, t_nd_gas_rel[0, 0]], color="tab:orange", zorder=1)
    axs[i].plot([s_gas_rel[0, n_lim], s_gas_rel[0, n_lim]], [t_nd_gas_rel[0, n_lim], t_3_gas_rel[0, n_lim]], color="tab:orange", zorder=1)

    axs[i].scatter(

        [s_liq_rel[0, 0], s_liq_rel[0, n_lim], 0, s_gas_rel[0, n_lim]],
        [t_nd_liq_rel[0, 0], t_nd_liq_rel[0, n_lim], t_nd_gas_rel[0, 0], t_3_gas_rel[0, n_lim]],
        s=30, facecolors='white', edgecolors=["tab:blue", "tab:blue", "tab:orange", "tab:orange"], zorder=2

    )

    axs[i].set_xlabel("$s\\ /\\ s_{max}$")
    axs[i].set_ylabel("$T\\ /\\ T_{max}$")
    axs[i].set_title(

        "${{ \\nabla T_{{rocks}} }}^{{\\#\\#}} = 10^{{ {:0.1f} }}$".format(

            np.log10(grad_nd_nd_list[n_lim])

        )

    )

    axs[i].legend()

    if i > 0:

        axs[i].get_yaxis().set_visible(False)



plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.show()