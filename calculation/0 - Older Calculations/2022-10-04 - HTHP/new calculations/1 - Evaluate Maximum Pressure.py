# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.ndimage import gaussian_filter, uniform_filter
from main_code.constants import CALCULATION_FOLDER
from UNISIMConnect import UNISIMConnector
from scipy.interpolate import griddata
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import numpy as np
import os


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "1 - Max Pressure Analysis.usc")
output_file = os.path.join(base_folder, "00 - Output", "2 - Test Different Temperature and Depth.csv")

depth_cell = "C3"
grad_cell = "C4"
force_ihx_cell = "C1"
p_steam_cell = "C8"
dt_max_cell = "C12"

p_max_cell = "C14"
ihx_on_cell = "C15"
t_ihx_out_cell = "D15"
ihx_power_cell = "C23"


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #
depth_list = np.round(np.linspace(600, 5000, 15), 1)
grad_list = np.round(np.linspace(15, 100, 15), 1)
p_steam_list = [

    {"ihx_off": 1, "P": 12},
    {"ihx_off": 1, "P": 3},
    {"ihx_off": 0, "P": 12}

]
X, Y = np.meshgrid(depth_list, grad_list, indexing='ij')

p_max_list = np.empty(np.append(X.shape, len(p_steam_list)))
ihx_on_list = np.empty(np.append(X.shape, len(p_steam_list)))
t_ihx_list = np.empty(np.append(X.shape, len(p_steam_list)))
ihx_power_list = np.empty(np.append(X.shape, len(p_steam_list)))

p_max_list[:] = np.nan
ihx_on_list[:] = np.nan
t_ihx_list[:] = np.nan
ihx_power_list[:] = np.nan

pbar = tqdm(desc="calculation", total=len(depth_list) * len(grad_list) * len(p_steam_list))
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("CALCULATION")

    for k in range(len(p_steam_list)):

        spreadsheet.set_cell_value(force_ihx_cell, p_steam_list[k]["ihx_off"])
        unisim.wait_solution()

        spreadsheet.set_cell_value(p_steam_cell, p_steam_list[k]["P"]*100)
        unisim.wait_solution()

        for i in range(X.shape[0]):

            spreadsheet.set_cell_value(depth_cell, X[i, 0])
            unisim.wait_solution()

            for j in range(X.shape[1]):

                spreadsheet.set_cell_value(grad_cell, Y[i, j])
                unisim.wait_solution()

                p_max_list[i, j, k] = spreadsheet.get_cell_value(p_max_cell)
                ihx_on_list[i, j, k] = spreadsheet.get_cell_value(ihx_on_cell)
                t_ihx_list[i, j, k] = spreadsheet.get_cell_value(t_ihx_out_cell)
                ihx_power_list[i, j, k] = spreadsheet.get_cell_value(ihx_power_cell)

                pbar.update(1)

pbar.close()
p_max_list = p_max_list / 100 # kPa to bar


# %%------------   REFINE RESULTS                         -----------------------------------------------------------> #
X_fine = np.linspace(np.min(X), np.max(X), 100)
Y_fine = np.linspace(np.min(Y), np.max(Y), 100)
X_fine, Y_fine = np.meshgrid(X_fine, Y_fine)
filter_gaussian = True

p_max_fine = np.empty(np.append(X_fine.shape, len(p_steam_list)))
ihx_on_fine = np.empty(np.append(X_fine.shape, len(p_steam_list)))
t_ihx_fine = np.empty(np.append(X_fine.shape, len(p_steam_list)))
ihx_power_fine = np.empty(np.append(X_fine.shape, len(p_steam_list)))

mask = ~np.isnan(p_max_list)

for k in range(len(p_steam_list)):

    points = np.vstack((X[mask[: , :, k]].ravel(), Y[mask[: , :, k]].ravel())).T
    p_max_clear = p_max_list[:, :, k]
    ihx_on_clear = ihx_on_list[:, :, k]
    t_ihx_clear = t_ihx_list[:, :, k]
    ihx_power_clear = ihx_power_list[:, :, k]

    p_max_fine[:, :, k] = griddata(points, p_max_clear[mask[: , :, k]], (X_fine, Y_fine), method='cubic')
    ihx_on_fine[:, :, k] = griddata(points, ihx_on_clear[mask[: , :, k]], (X_fine, Y_fine), method='cubic')
    t_ihx_fine[:, :, k] = griddata(points, t_ihx_clear[mask[: , :, k]], (X_fine, Y_fine), method='cubic')
    ihx_power_fine[:, :, k] = griddata(points, ihx_power_clear[mask[: , :, k]], (X_fine, Y_fine), method='cubic')

    if filter_gaussian:
        sigma_val = 3
        p_max_fine[:, :, k] = gaussian_filter(p_max_fine[:, :, k], sigma=sigma_val)
        ihx_on_fine[:, :, k] = gaussian_filter(ihx_on_fine[:, :, k], sigma=sigma_val)
        t_ihx_fine[:, :, k] = gaussian_filter(t_ihx_fine[:, :, k], sigma=sigma_val)
        ihx_power_fine[:, :, k] = gaussian_filter(ihx_power_fine[:, :, k], sigma=sigma_val)

    else:
        p_max_fine[:, :, k] = uniform_filter(p_max_fine[:, :, k], size=5)
        ihx_on_fine[:, :, k] = uniform_filter(ihx_on_fine[:, :, k], size=5)
        t_ihx_fine[:, :, k] = uniform_filter(t_ihx_fine[:, :, k], size=5)
        ihx_power_fine[:, :, k] = uniform_filter(ihx_power_fine[:, :, k], size=5)

    ihx_power_fine[(ihx_power_fine[:, :, k] < np.max(ihx_power_fine[:, :, k]) * 0.05), k] = 0


# %%------------   PLOT CONTOUR AND IHX ON                -----------------------------------------------------------> #
# <------------------ INPUT DEFINITION -------------------------------------------->
# k_list = np.array([[0, 1, 2]])
# plot_ihx_power = np.array([[False, False, False]])
# hide_colorbar = True
# size_fact = {"col":5, "row":4.5}

k_list = np.array([[-1, 2]])
plot_ihx_power = np.array([[True, True]])
hide_colorbar = False
size_fact = {"col":6.5, "row":4.5}

tight_pad = 2
label_font_size = 11
title_font_size = 12
contour_label_font_size = 10
contour_label_color = 'black'

n_levels = 11
spec_levels_ihx_on = [100, 125, 150, 175, 200, 225, 250]
spec_levels_ihx_off = [10, 50, 100, 250, 500, 1000]
spec_levels_all = [10, 50, 100, 200, 1000]

background_alpha = 0.5
background_cmap = 'viridis'


# <------------------ CALCULATION      -------------------------------------------->
n_rows = k_list.shape[0]
n_cols = k_list.shape[1]
fig, axes = plt.subplots(

    nrows=n_rows, ncols=n_cols,
    figsize=(n_cols*size_fact["col"], n_rows*size_fact["row"])

)

axes = np.array([axes])

k_list = k_list.flatten()
axes = axes.flatten()
plot_ihx_power = plot_ihx_power.flatten()

for i in range(n_rows * n_cols):

    k = k_list[i]

    if k >= 0:

        if p_steam_list[k]["ihx_off"] == 0:

            ihx_str = "- with\ IHX "
            specified_levels = spec_levels_ihx_on

        else:

            ihx_str = ""
            specified_levels = spec_levels_ihx_off

        contour_plot = axes[i].contour(

            X_fine, Y_fine, p_max_fine[:,:,k],
            levels=specified_levels, colors=contour_label_color

        )
        axes[i].clabel(

            contour_plot, inline=True,
            fontsize=contour_label_font_size,
            fmt=lambda x: f'{x:.0f} bar'

        )

        levels = np.linspace(0, 1, n_levels)
        if p_steam_list[k]["ihx_off"] == 0 and plot_ihx_power[i]:

            z_rel = ihx_power_fine / np.max(ihx_power_fine[:,:,k])
            cbar_label = "IHX power (rel) [%]"

            cont = axes[i].contourf(

                X_fine, Y_fine, z_rel[:,:,k],
                levels=levels, alpha=background_alpha,
                extend='both', cmap=background_cmap + '_r'
            )

            where = z_rel[:, :, k] == 0
            axes[i].contourf(X_fine, Y_fine, where, levels=[0.5, 1], colors="white")

        else:

            z_rel = p_max_fine[:, :, k] / np.max(p_max_fine[:, :, k])
            cbar_label = "Pressure [%]"

            cont = axes[i].contourf(

                X_fine, Y_fine, z_rel, levels=levels,
                cmap=background_cmap, extend='both',
                alpha=background_alpha

            )

        if not hide_colorbar:

            cbar = fig.colorbar(

                cont, ax=axes[i],
                orientation='vertical',
                fraction=0.046, pad=0.04,

            )
            cbar.ax.tick_params(labelsize=label_font_size)
            cbar.set_label(cbar_label, fontsize=label_font_size)

        axes[i].set_title(

            "$\it{{(Steam @ {:.0f} bar){}}}$".format(

                p_steam_list[k]["P"], ihx_str

            ), fontsize=title_font_size, loc='center'

        )

    else:

        specified_levels = spec_levels_all

        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        colors = cycle(default_colors)
        colors_without = cycle(default_colors)

        legend_entries = list()
        legend_labels = list()
        ovr_base_label = 'Steam @ {:.0f} bar\n{}'
        for n in range(len(p_steam_list)):

            if p_steam_list[n]["ihx_off"] == 0:

                ihx_str = "$\it{with\ IHX}$"
                linestyle_curr = "-"
                color = next(colors)

            else:

                ihx_str = ("$\it{without\ IHX}$")
                linestyle_curr = "--"
                color = next(colors_without)

            contour_plot = axes[i].contour(

                X_fine, Y_fine, p_max_fine[:, :, n],
                levels=specified_levels, colors=color,
                linestyles=linestyle_curr

            )

            axes[i].clabel(contour_plot, inline=True, fontsize=contour_label_font_size, fmt=lambda x: f'{x:.0f} bar')
            legend_entries.append(Line2D([0], [0], color=color, lw=1, linestyle=linestyle_curr))
            legend_labels.append(ovr_base_label.format(p_steam_list[n]["P"], ihx_str))

        axes[i].legend(

            legend_entries, legend_labels,
            fontsize=label_font_size, loc='center left', bbox_to_anchor=(1, 0.5),

        )

        # axes[i].set_ylim([np.min(X_fine), np.max(Y_fine)])
        # axes[i].set_xlim([np.min(X_fine), np.max(Y_fine)])

    axes[i].set_xlabel('depth [m]', fontsize=label_font_size, loc='right')
    axes[i].set_ylabel('grad [°C/km]', fontsize=label_font_size, loc='top')
    axes[i].tick_params(axis='both', which='major', labelsize=label_font_size)

plt.tight_layout(pad = tight_pad)
plt.show()