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
a_curr = 15
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)
dz_label = "${{\\Delta z}}^{{\\#}} = 10^{{ {:0.1f} }}$"
line_styles = ["-", "--", "-.", "-"]
w_dot_rel = ovr_res["w_dot_nds"] / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

lines = list()
for i in range(n_t_rel):

    lines = list()
    for k in range(n_depth):

        lines.append(

            ax.plot(

                grad_nd_nd[:, k], ovr_res["w_dot_nds"][a_curr, i, :, k], #line_styles[i],
                label=dz_label.format(np.log10(dz_nd_list[k])),
                color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

            )[0]

        )

        ax.plot(

            grad_nd_nd[:, k], spc_work_liq[:, k], #line_styles[i],
            label=dz_label.format(np.log10(dz_nd_list[k])),
            color=cmap(norm((k + 1) / (len(dz_nd_list) + 1)))

        )

        ax.plot(

            grad_nd_nd[:, k], spc_work_gas[:, k], #line_styles[i],
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
a_curr = 15
fig, axs = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(13, 5)
cmap = plt.get_cmap('viridis')
t_rel_label = "$T_{{rel}} = 10^{{ {:0.1f} }}$"

w_dot_rel = ovr_res["w_dot_nds"] / spc_work_liq
dw_dot_rel = w_dot_rel[:, 1:, :] - w_dot_rel[:, :-1, :]

lines = list()
min_grad_dz_all = list()
max_grad_dz_all = list()
max_grad_dz_ovr_all = list()

t_rel_list = [0.5, 0.75, 1.1, 2]

sub_t_rel_list = t_rel_list[:-1]
for i in range(len(sub_t_rel_list)):

    min_grad = list()
    max_grad = list()
    max_grad_ovr = list()

    min_grad_dz = list()
    max_grad_dz = list()
    max_grad_dz_ovr = list()

    for k in range(len(dz_nd_list)):

        try:

            max_dw_dot_rel = np.nanmax(dw_dot_rel[a_curr, i, :, k])
            max_w_dot_rel = np.nanmax(w_dot_rel[a_curr, i, np.where(dw_dot_rel[a_curr, i, :, k] > 0), k])
            min_w_dot_rel = np.nanmin(w_dot_rel[a_curr, i, np.where(dw_dot_rel[a_curr, i, :, k] > 0.1 * max_dw_dot_rel), k])

            min_j = np.where(w_dot_rel[a_curr, i, :, k] == min_w_dot_rel)[0][0]
            max_j = np.where(dw_dot_rel[a_curr, i, :, k] == max_dw_dot_rel)[0][0]
            max_ovr_j = np.where(w_dot_rel[a_curr, i, :, k] == max_w_dot_rel)[0][0]

            min_grad.append(grad_nd_nd[min_j, k])
            max_grad.append(grad_nd_nd[max_j, k])
            max_grad_ovr.append(grad_nd_nd[max_ovr_j, k])

            min_grad_dz.append(grad_nd_nd[min_j, k] * dz_nd[min_j, k])
            max_grad_dz.append(grad_nd_nd[max_j, k] * dz_nd[max_j, k])
            max_grad_dz_ovr.append(grad_nd_nd[max_ovr_j, k] * dz_nd[max_ovr_j, k])

        except:

            min_grad.append(np.nan)
            max_grad.append(np.nan)
            max_grad_ovr.append(np.nan)

            min_grad_dz.append(np.nan)
            max_grad_dz.append(np.nan)
            max_grad_dz_ovr.append(np.nan)

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


# %%-------------------------------------   PLOT CONTOUR                        -------------------------------------> #
fig, ax = plt.subplots(

    1, 1, figsize=(12, 5),
    dpi=300

)

n_level = 30
n_isoline_level = 6
isoline_width = 0.75
isoline_alpha = 0.6
plot_contour_lines = False

max_log = 2
min_log = np.log10(np.min(res_tp_mesh))

max_log_label = 2
min_log_label = -1

levels = np.logspace(min_log, max_log, n_level)
ticks_levels = np.logspace(min_log_label, max_log_label, (max_log_label - min_log_label) + 1)
isoline_levels = np.logspace(min_log,max_log_label, n_isoline_level)

cs = ax.contourf(

    v_tv_mesh, t_tv_mesh,
    res_tv_mesh, levels=levels,
    zorder=1, norm=colors.LogNorm()

)

thinks_labels = list()
for ticks_level in ticks_levels:

    found = False
    for i in range(3):

        if ticks_level >= 10**-i:
            thinks_labels.append("{{:.{}f}}".format(i).format(ticks_level))
            found = True
            break

    if not found:
        thinks_labels.append("$10^{{ {:.0f} }}$".format(np.log10(ticks_level)))

axs[2].get_xaxis().set_visible(False)
cbar = fig.colorbar(cs, cax=axs[2])
cbar.set_ticks(ticks_levels)
cbar.set_ticklabels(thinks_labels)
cbar.ax.set_ylabel("${\\nabla T_{rocks}}^{\\#}_{lim}$ [-]")

if plot_contour_lines:

    cl_p = axs[0].contour(

        p_tp_mesh, t_tp_mesh, res_tp_mesh,
        levels=isoline_levels, colors="black",
        alpha=isoline_alpha, linewidths=isoline_width,
        zorder=1

    )
    cl_v = axs[1].contour(

        v_tv_mesh, t_tv_mesh, res_tv_mesh,
        levels=isoline_levels, colors="black",
        alpha=isoline_alpha, linewidths=isoline_width,
        zorder=1

    )
    labels_p = axs[0].clabel(cl_p, inline=0.1, fontsize=7, zorder=1)
    labels_v = axs[1].clabel(cl_v, inline=0.1, fontsize=7, zorder=1)

    for l in labels_p + labels_v:
        l.set_rotation(0)

axs[0].plot(p_sat_arr[:, 0], t_rels, color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:-1, 0], t_rels[1:], color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:, 1], t_rels[:], color="black", linewidth=2, zorder=2)

plots_each_max = 20
tmp_list = np.linspace(0, len(v_max[0]) - 1, len(v_max[0]))
mark_indices = np.where(np.mod(tmp_list, plots_each_max) == 0)
axs[0].plot(p_max[1], p_max[0], "-", color="tab:orange", zorder=3)
axs[1].plot(v_max[1], v_max[0], "-", color="tab:orange", zorder=3)

axs[0].scatter(

    p_max[1, mark_indices], p_max[0, mark_indices],
    s=20, edgecolors="tab:orange",
    facecolors='white', zorder=100

)
axs[1].scatter(

    v_max[1, mark_indices], v_max[0, mark_indices],
    s=20, edgecolors="tab:orange",
    facecolors='white', zorder=100

)


for ax in axs[:-1]:

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("$T_{rel}$ [-]")

axs[0].set_xlabel("$p_{rel}$ [-]")
axs[1].set_xlabel("$(v - b) / (v_{crit} - b)$ [-]")
axs[0].set_xlim((np.min(p_rels), np.max(p_rels)))
axs[1].set_xlim((np.min(v_rels), np.max(v_rels)))

plt.tight_layout(pad=1)
plt.show()
