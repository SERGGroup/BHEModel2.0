# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from REFPROPConnector import ThermodynamicPoint
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   CALCULATIONS OPTIONS                -------------------------------------> #
n_points = 75

t_rels = np.logspace(np.log10(0.2), np.log10(2), n_points)[:-1]
v_rels = np.logspace(-2, 3, n_points)
p_rels = np.logspace(-2, 2, n_points)

t_tv_mesh, v_tv_mesh = np.meshgrid(t_rels, v_rels, indexing='ij')
t_tp_mesh, p_tp_mesh = np.meshgrid(t_rels, p_rels, indexing='ij')


# %%-------------------------------------   INIT FLUID                          -------------------------------------> #
fluid = "Methane"
tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
bhe_in = ThermodynamicPoint([fluid], [1])
tp.set_variable("T", tp.RPHandler.TC)
tp.set_variable("P", tp.RPHandler.PC)


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
res_tp_mesh = np.empty(p_tp_mesh.shape)
res_tv_mesh = np.empty(v_tv_mesh.shape)
res_tp_mesh[:] = np.nan
res_tv_mesh[:] = np.nan
p_sat_arr = np.zeros((len(t_rels), 1))
v_sat_arr = np.zeros((len(t_rels), 2))

pbar = tqdm(desc="Calculating Points", total=n_points * len(t_rels) * 2)

for i in range(len(t_rels)):

    t_curr = crit_state.t * t_rels[i]
    sat_states = fluid_curr.get_sat_state(t=t_curr)
    p_sat_arr[i, 0] = sat_states[0].p
    v_sat_arr[i, 0] = sat_states[0].v
    v_sat_arr[i, 1] = sat_states[1].v

    for j in range(len(v_rels)):

        v_curr = v_rels[j] * (crit_state.v - fluid_curr.b) + fluid_curr.b
        curr_state = fluid_curr.get_state(v=v_curr, t=t_curr, sat_states=sat_states)

        if not curr_state.bifase:
            res_tv_mesh[i, j] = -curr_state.t / curr_state.v * curr_state.dpdt / curr_state.dpdv

        pbar.update(1)

    for k in range(len(p_rels)):

        p_curr = p_rels[k] * crit_state.p
        curr_state = fluid_curr.get_state(p=p_curr, t=t_curr, sat_states=sat_states)

        if not curr_state.bifase:
            res_tp_mesh[i, k] = -curr_state.t / curr_state.v * curr_state.dpdt / curr_state.dpdv

        pbar.update(1)


p_sat_arr = p_sat_arr / fluid_curr.p_crit
v_sat_arr = (v_sat_arr - fluid_curr.b) / (fluid_curr.v_crit - fluid_curr.b)

pbar.close()

# %%-------------------------------------   IDENTIFY MAX VALUES                 -------------------------------------> #
v_max = [list(), list()]
p_max = [list(), list()]
firs_j = False
firs_k = False
for i in range(len(t_rels)):

    if t_rels[i] > 1:

        max_r_v_values = np.max(res_tv_mesh[i, :])
        j = np.where(res_tv_mesh[i, :] == max_r_v_values)

        if not firs_j:

            v_max[0].append(t_tv_mesh[i, j])
            v_max[1].append(v_tv_mesh[i, j])

            if j[0][0] == len(v_rels) - 1 or j[0][0] == 0:

                firs_j = True

        max_r_t_values = np.max(res_tp_mesh[i, :])
        k = np.where(res_tp_mesh[i, :] == max_r_t_values)

        if not firs_k:

            p_max[0].append(t_tp_mesh[i, k])
            p_max[1].append(p_tp_mesh[i, k])

            if k[0][0] == len(p_rels) - 1 or k[0][0] == 0:

                firs_k = True

v_max = np.array(v_max)[:, :, 0, 0]
p_max = np.array(p_max)[:, :, 0, 0]


# %%-------------------------------------   PLOT CONTOUR                        -------------------------------------> #
c_bar_ratio = 15
fig, axs = plt.subplots(

    1, 3, figsize=(12, 5),
    dpi=300, gridspec_kw={'width_ratios': [c_bar_ratio, c_bar_ratio, 1]}

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

axs[0].contourf(

    p_tp_mesh, t_tp_mesh,
    res_tp_mesh, levels=levels,
    zorder=1, norm=colors.LogNorm()

)

cs = axs[1].contourf(

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


# %%-------------------------------------   PLOT REGIONS                        -------------------------------------> #
fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

isoline_width = 1
isoline_alpha = 1
plot_contour_lines = True
levels = [np.min(res_tp_mesh), 0.98, 1.2, np.max(res_tp_mesh)]
isoline_levels = [0.98, 1.2]
cmap = mpl.colormaps["cividis"]
norm = colors.Normalize(vmin=0.98, vmax=1.2)

axs[0].contourf(

    p_tp_mesh, t_tp_mesh,
    res_tp_mesh, levels=levels,
    zorder=1, norm=norm,
    cmap=cmap

)

cs = axs[1].contourf(

    v_tv_mesh, t_tv_mesh,
    res_tv_mesh, levels=levels,
    zorder=1, norm=norm,
    cmap=cmap

)

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

axs[0].plot(p_sat_arr[:, 0], t_rels, color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:-1, 0], t_rels[1:], color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:-1, 1], t_rels[1:], color="black", linewidth=2, zorder=2)

for ax in axs:

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("$T_{rel}$ [-]")

axs[0].set_xlabel("$p_{rel}$ [-]")
axs[1].set_xlabel("$(v - b) / (v_{crit} - b)$ [-]")
axs[0].set_xlim((np.min(p_rels), np.max(p_rels)))
axs[1].set_xlim((np.min(v_rels), np.max(v_rels)))

props = dict(boxstyle='round', facecolor='white', alpha=1)
axs[0].text(0.35, 0.75, "A", transform=axs[0].transAxes, fontsize=14, bbox=props)
axs[0].text(0.75, 0.15, "B", transform=axs[0].transAxes, fontsize=14, bbox=props)
axs[0].text(0.6, 0.5, "C", transform=axs[0].transAxes, fontsize=14, bbox=props)

axs[1].text(0.75, 0.75, "A", transform=axs[1].transAxes, fontsize=14, bbox=props)
axs[1].text(0.15, 0.55, "B", transform=axs[1].transAxes, fontsize=14, bbox=props)
axs[1].text(0.45, 0.5, "C", transform=axs[1].transAxes, fontsize=14, bbox=props)

plt.tight_layout(pad=1)
plt.show()

