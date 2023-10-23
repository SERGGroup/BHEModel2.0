# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   INPUT POINTS                        -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "0.support", "res"

)


def init_result_dict(res_shape):

    result_dict = {

        "grad_nd_rocks": np.empty(res_shape),
        "grad_rocks": np.empty(res_shape),
        "depth": np.empty(res_shape),
        "alpha_in": np.empty(res_shape),
        "t_rocks": np.empty(res_shape),
        "w_dot_nds": np.empty(res_shape),
        "ex_dot_nds": np.empty(res_shape),
        "eta_exs": np.empty(res_shape),
        "w_dex_mins": np.empty(res_shape),
        "w_dex_maxs": np.empty(res_shape)

    }

    for key in result_dict.keys():
        result_dict[key][:] = np.nan

    return result_dict


test_arr = np.load(os.path.join(result_folder, "w_dot_nds.npy"), allow_pickle=True)
ovr_res = init_result_dict(test_arr.shape)
for key in ovr_res.keys():
    ovr_res[key] = np.load(os.path.join(result_folder, "{}.npy".format(key)), allow_pickle=True)

filename = os.path.join(result_folder, "elapsed_times.npy")
times = np.load(filename, allow_pickle=True)

n_p_rel = test_arr.shape[0]
n_t_rel = test_arr.shape[1]
n_grad = test_arr.shape[2]
n_depth = test_arr.shape[3]


# %%-------------------------------------   EVALUATE LIQUID AND IDEAL GASSES    -------------------------------------> #
dz_nd_list = np.logspace(-2, 0, n_depth)
grad_nd_nd_list = np.logspace(-2, 4.5, n_grad)
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


# %%-------------------------------------   IDENTIFY INITIAL RISE               -------------------------------------> #
a_curr = 2
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)
dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
line_styles = ["-", "--", "-.", "-"]
w_dot_rel = ovr_res["w_dot_nds"] / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

for i in range(n_t_rel):

    lines = list()
    for k in range(n_depth):

        lines.append(

            ax.plot(

                grad_nd_nd[:, k], ovr_res["w_dot_nds"][a_curr, i, :, k], line_styles[i],
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

            grad_nd_nd[:, k], spc_work_gas[:, k], line_styles[i],
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
