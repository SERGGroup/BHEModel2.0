# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from matplotlib.colors import Normalize, LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, uniform_filter
from main_code.constants import CALCULATION_FOLDER
from matplotlib.ticker import FormatStrFormatter
from UNISIMConnect import UNISIMConnector
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tqdm import tqdm
import pandas as pd
import numpy as np
import os


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "3 - Optimization-better.usc")
output_file = os.path.join(base_folder, "00 - Output", "2 - Optimization.csv")

t_sg_perc_cell = "C12"
sep_perc_cell = "C13"

m_ratio_cell = "G4"
m_ratio_rel_cell = "G6"
COP_cell = "C17"
eta_el_cell = "C18"
eta_cell = "C19"

w_net_cell = "G13"
q_steam_cell = "G17"

t_set_cell = "E14"
pinch_error_cell = "E15"

rel_power_cells = {

    "main_turb":     "H9",
    "hp_turb":      "H10",
    "comp":         "H11",
    "pump":         "H12",
    "q_steam":      "H17"

}

m_ratio_cells = {

    "m_ratio":      "H4",
    "m_ratio_rel":  "H6",

}

def evaluate_params(sheet, sep_perc: np.ndarray, use_rel_ratio: bool = True):

    net_power = sep_perc * sheet.get_cell_value(rel_power_cells["main_turb"])
    net_power += (1 - sep_perc) * sheet.get_cell_value(rel_power_cells["hp_turb"])
    net_power -= (1 - sep_perc) * sheet.get_cell_value(rel_power_cells["comp"])
    net_power -= (1 - sep_perc) * sheet.get_cell_value(rel_power_cells["pump"])
    net_power = net_power / ((1 - sep_perc) * sheet.get_cell_value(rel_power_cells["q_steam"]))

    if use_rel_ratio:

        m_ratio = (1 - sep_perc) * sheet.get_cell_value(m_ratio_cells["m_ratio_rel"])

    else:

        m_ratio = (1 - sep_perc) * sheet.get_cell_value(m_ratio_cells["m_ratio"])

    return net_power, m_ratio


results = {

    t_sg_perc_cell: list(),
    sep_perc_cell: list(),

    m_ratio_cell: list(),
    m_ratio_rel_cell: list(),
    COP_cell: list(),
    eta_el_cell: list(),
    eta_cell: list(),

}

for i in range(9, 14, 1):

    results.update({"G{}".format(i): list()})

for i in range(16, 20, 1):

    results.update({"G{}".format(i): list()})


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #
t_sg_perc_list = np.linspace(0, 1, 45)
t_sg_perc_list[0] = 0.001

sep_perc_list = np.linspace(0, 1, 45)
sep_perc_list[0] = 0.01
sep_perc_list[-1] = 1 - 0.01

X, Y = np.meshgrid(t_sg_perc_list, sep_perc_list, indexing='ij')

COP_list = np.empty(X.shape)
eta_list = np.empty(X.shape)
m_ratio_list = np.empty(X.shape)
w_rel_list = np.empty(X.shape)
q_steam_list = np.empty(X.shape)

COP_list[:] = np.nan
eta_list[:] = np.nan
m_ratio_list[:] = np.nan
w_rel_list[:] = np.nan
q_steam_list[:] = np.nan

for key in results.keys():
    results[key] = list()

pbar = tqdm(desc="calculation", total=len(t_sg_perc_list))
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("CALCULATION")

    for i in range(X.shape[0]):

        spreadsheet.set_cell_value(t_sg_perc_cell, X[i, 0])
        unisim.wait_solution()

        t_ihx_per_bound = np.array([0.2, 1])
        for j in range(30):

            t_ihx_per = np.mean(t_ihx_per_bound)
            spreadsheet.set_cell_value(t_set_cell, t_ihx_per)
            unisim.wait_solution()

            error = spreadsheet.get_cell_value(pinch_error_cell)

            try:

                if abs(error) < 0.001:

                    break

            except:

                t_ihx_per_bound[0] = t_ihx_per

            else:

                if error > 0:
                    t_ihx_per_bound[0] = t_ihx_per

                else:
                    t_ihx_per_bound[1] = t_ihx_per

        w_rel_list[i, :], m_ratio_list[i, :] = evaluate_params(sheet=spreadsheet, sep_perc=sep_perc_list, use_rel_ratio=False)
        pbar.update(1)

pbar.close()
pd.DataFrame(results).to_csv(output_file)


# %%------------   REFINE RESULTS and IDENTIFY PARETO     -----------------------------------------------------------> #
X_fine = np.linspace(np.min(X), np.max(X), 150)
Y_fine = np.linspace(np.min(Y), np.max(Y), 150)
# X_fine = np.linspace(0.2, 0.8, 150)
# Y_fine = np.linspace(0.2, 0.8, 150)
X_fine, Y_fine = np.meshgrid(X_fine, Y_fine)

filter_gaussian = True
sigma_val = 1
method = 'cubic'
points = np.vstack((X.ravel(), Y.ravel())).T

COP_fine = griddata(points, COP_list.ravel(), (X_fine, Y_fine), method=method)
eta_fine = griddata(points, eta_list.ravel(), (X_fine, Y_fine), method=method)
m_ratio_fine = griddata(points, m_ratio_list.ravel(), (X_fine, Y_fine), method=method)
q_steam_fine = griddata(points, q_steam_list.ravel(), (X_fine, Y_fine), method=method)
w_rel_fine = griddata(points, w_rel_list.ravel(), (X_fine, Y_fine), method=method)

if filter_gaussian:

    COP_fine = gaussian_filter(COP_fine, sigma=sigma_val)
    eta_fine = gaussian_filter(eta_fine, sigma=sigma_val)
    m_ratio_fine = gaussian_filter(m_ratio_fine, sigma=sigma_val)
    q_steam_fine = gaussian_filter(q_steam_fine, sigma=sigma_val)
    w_rel_fine = gaussian_filter(w_rel_fine, sigma=sigma_val)

pareto_points = {

    "X": list(),
    "Y": list(),
    "Z": list(),

    "m_ratio": list(),
    "w_net": list()

}
m_ratio_flat = m_ratio_fine.flatten()
w_rel_flat = w_rel_fine.flatten()
t_sg_perc_flat = X_fine.flatten()
sep_perc_flat = Y_fine.flatten()
otm_points = np.vstack((m_ratio_flat, w_rel_flat)).T

pbar = tqdm(desc="Searching for Pareto Front Points", total=len(m_ratio_flat))
for i, curr_point in enumerate(otm_points):

    is_pareto = True
    for j, other_point in enumerate(otm_points):

        if all(other_point >= curr_point) and any(other_point > curr_point):
            is_pareto = False
            break

    pbar.update(1)
    if is_pareto:
        Z = 1 / (1 + np.exp(((w_rel_flat[i] / m_ratio_flat[i])) / 20))
        pareto_points["X"].append(t_sg_perc_flat[i])
        pareto_points["Y"].append(sep_perc_flat[i])
        pareto_points["Z"].append(Z)
        pareto_points["m_ratio"].append(m_ratio_flat[i])
        pareto_points["w_net"].append(w_rel_flat[i])

pbar.close()


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
# <------------------ INPUT DEFINITION   ------------------------------------------>
a = 1
omega = 0.55
sigma_val = 5
third_levels = 30

window_size = 15  # Should be an odd number
poly_order = 2   # Polynomial order

x_min_range = np.log(w_rel_list + 1) - omega * np.log(1 / m_ratio_list)

title_fontsize = 14
label_font_size = 10
tight_pad = 1
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.5))
x_label = '$T_{sg\ \%}$ [-]'
y_label = '$sep_{\%}$ [-]'

# <------------------ FIRST AX          ------------------------------------------>
axes[0].set_title('$W_{net\ rel}$', fontsize=title_fontsize)

negative_levels = np.linspace(-0.5, 0, 50)
positive_levels = np.linspace(0, 5, 50)
combined_levels = np.concatenate([negative_levels, positive_levels[1:]])

# Create a custom colormap that blends two colormaps
# ('viridis' for negative values and 'plasma_r' for positive)
cmap_neg = cm.get_cmap('viridis', len(negative_levels))
cmap_pos = cm.get_cmap('plasma_r', len(positive_levels))

cmap_neg = cmap_neg(np.linspace(0, 1, len(negative_levels)))
cmap_pos = cmap_pos(np.linspace(0, 1, len(negative_levels)))

min_value = np.min(negative_levels)
max_value = np.max(positive_levels)
d_value = max_value - min_value
norm = Normalize(vmin=min_value, vmax=max_value)

colors_combined = list()
for i in range(len(negative_levels)):
    colors_combined.append(( (negative_levels[i]-min_value)/d_value, cmap_neg[i]))

for j in range(len(positive_levels) - 1):
    colors_combined.append(( (positive_levels[j+1]-min_value)/d_value, cmap_pos[j+1]))

# Combine the two colormaps into one
combined_cmap = LinearSegmentedColormap.from_list('combined_cmap', colors_combined)

contour_plot = axes[0].contour(X_fine, Y_fine, w_rel_fine, levels=[0], colors='black', linewidths=2)
contourf_plot = axes[0].contourf(X_fine, Y_fine, w_rel_fine, levels=combined_levels, cmap=combined_cmap, norm=norm, extend='both')
axes[0].clabel(

    contour_plot, inline=True,
    fontsize=10,
    fmt=lambda x: f'$W_{{net\ rel}}$ = {x:.0f}'

)

cbar = fig.colorbar(contourf_plot, label='$W_{net\ rel}$ [-]', ax=axes[0])
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


# <------------------ SECOND AX          ------------------------------------------>
axes[1].set_title('$m_{ratio}$', fontsize=title_fontsize)
contourf_plot = axes[1].contourf(X_fine, Y_fine, m_ratio_fine, levels=25, cmap="viridis")
cbar = fig.colorbar(contourf_plot, label='$m_{ratio}$ [-]', ax=axes[1], extend='both')
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# <------------------ FINAL CALCULATIONS ------------------------------------------>
for ax in axes:

    ax.set_xlabel(x_label, fontsize=label_font_size, loc='right')
    ax.set_ylabel(y_label, fontsize=label_font_size, loc='top')
    ax.tick_params(axis='both', which='major', labelsize=label_font_size)

plt.tight_layout(pad = tight_pad)
plt.show()


# %%------------   PLOT PARETO                           -----------------------------------------------------------> #
# <------------------ INPUT DEFINITION   ------------------------------------------>
sigma_val = 2
window_size = 5
omega = 0.45
window_size_filter = 5  # Should be an odd number
poly_order = 2   # Polynomial order

gaussian_filter_before = True
gaussian_filter_after = False
use_savgol = True
c_smooth = "black"
c_curve = "tab:orange"
c_pareto = "tab:blue"

omega_list = np.linspace(0, 1, 150)
omega_values = np.linspace(0, 1, 21)[1:-1]
levels_minx = np.linspace(0, 0.5, 25)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.5))


# <------------------ FIRST AX           ------------------------------------------>
axes[0].set_title('$x_{{min}}$ - $\Omega$={}'.format(omega), fontsize=title_fontsize)
a = np.log(w_rel_fine + 1)
b = np.log(1 / m_ratio_fine)

if gaussian_filter_before:
    a = gaussian_filter(a, sigma=sigma_val)
    b = gaussian_filter(b, sigma=sigma_val)

dx_smooth = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=1, delta=1)
dy_smooth = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=0, delta=1)
dtot_smoot_1 = dx_smooth * dx_smooth + dy_smooth * dy_smooth

if gaussian_filter_after:
    dtot_smoot_1 = gaussian_filter(dtot_smoot_1, sigma=sigma_val)

min_der = np.nanmin(dtot_smoot_1)
i_min = np.where(dtot_smoot_1 == min_der)
x_star = np.nanmean(X_fine[i_min])
y_star = np.nanmean(Y_fine[i_min])

z_value = (1 - omega) * a - omega * b
# z_value = dtot_smoot_1
z_value = (z_value - np.nanmin(z_value)) / (np.nanmax(z_value) - np.nanmin(z_value))

contourf_plot = axes[0].contourf(X_fine, Y_fine, z_value, levels=levels_minx, cmap="viridis", extend='both')
cbar = fig.colorbar(contourf_plot, label='$x_{min}$ [-]', ax=axes[0])
axes[0].plot(x_star, y_star, marker='*', markersize=15, color='#FFD700')

# <------------------ SECOND AX         ------------------------------------------>
axes[1].set_title('Pareto Front Vs. $\Omega$ Optimization'.format(omega), fontsize=title_fontsize)
a = np.log(w_rel_fine + 1)
b = np.log(1 / m_ratio_fine)

x_min_list = np.empty(omega_list.shape)
y_min_list = np.empty(omega_list.shape)
min_der_list = np.empty(omega_list.shape)
x_min_list[:] = np.nan
y_min_list[:] = np.nan
min_der_list[:] = np.nan

for i, omega in enumerate(omega_list):

    if use_savgol:

        dx = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=1, delta=1)
        dy = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=0, delta=1)
        x_curr = X_fine
        y_curr = Y_fine

    else:

        da_dx = a[1:, 1:] - a[:-1, 1:]
        da_dy = a[1:, 1:] - a[1:, :-1]
        db_dx = b[1:, 1:] - b[:-1, 1:]
        db_dy = b[1:, 1:] - b[1:, :-1]

        if gaussian_filter_before:
            da_dx = gaussian_filter(da_dx, sigma=sigma_val)
            da_dy = gaussian_filter(da_dy, sigma=sigma_val)
            db_dx = gaussian_filter(db_dx, sigma=sigma_val)
            db_dy = gaussian_filter(db_dy, sigma=sigma_val)

        dx = (1 - omega) * da_dx - omega * db_dx
        dy = (1 - omega) * da_dy - omega * db_dy
        x_curr = X_fine[1:, 1:]
        y_curr = Y_fine[1:, 1:]

    dtot = dx * dx + dy * dy

    if gaussian_filter_after:

        dtot = gaussian_filter(dtot, sigma=sigma_val)

    min_der = np.nanmin(dtot)
    i_min = np.where(dtot == min_der)

    try:

        x_min_list[i] = np.nanmean(x_curr[i_min])
        y_min_list[i] = np.nanmean(y_curr[i_min])
        min_der_list[i] = min_der

    except:

        print(i_min, i)

# Delete not minima
min_der_list = (min_der_list - np.nanmin(min_der_list)) / (np.nanmax(min_der_list) - np.nanmin(min_der_list))
not_minima = min_der_list > 4e-3

x_min_list[not_minima] = np.nan
y_min_list[not_minima] = np.nan

# Find the indices for omega = 0.25, 0.5, and 0.75
indices = [np.argmin(np.abs(omega_list - omega)) for omega in omega_values]

nan_mask = np.isnan(y_min_list)

x_min_series = pd.Series(x_min_list)
x_min_filled = x_min_series.fillna(method='ffill')
x_min_filled = x_min_filled.fillna(method='bfill')
x_min_cleaned = x_min_filled.to_numpy()
x_smooth = np.convolve(x_min_cleaned, np.ones(window_size) / window_size, mode='same')

y_min_series = pd.Series(y_min_list)
y_min_filled = y_min_series.fillna(method='ffill')
y_min_filled = y_min_filled.fillna(method='bfill')
y_min_cleaned = y_min_filled.to_numpy()
y_smooth = np.convolve(y_min_cleaned, np.ones(window_size) / window_size, mode='same')

# Plot Pareto points and the curve
axes[1].scatter(pareto_points["X"], pareto_points["Y"], c=c_pareto, alpha=0.1)
axes[1].plot(x_min_list, y_min_list, c=c_curve, linewidth=2, alpha=0.75)
axes[1].plot(x_smooth[~nan_mask], y_smooth[~nan_mask], c=c_smooth, linewidth=5, alpha=0.75)

# Calculate the slope (dx, dy) between consecutive points to determine the curve's direction
dx_list = np.diff(x_smooth)
dy_list = np.diff(y_smooth)

# Normalize the slope vectors to unit length (dx, dy)
norms = np.sqrt(dx_list**2 + dy_list**2)
dx_unit = dx_list / norms
dy_unit = dy_list / norms

# Add small perpendicular lines at specific omega values
for idx, omega in zip(indices, omega_values):

    # Find the normalized perpendicular direction to the curve at this point
    perp_dx = -dy_unit[idx-1]  # Perpendicular to the curve direction (rotate by 90 degrees)
    perp_dy = dx_unit[idx-1]

    # Define a small line length for the tick
    tick_length = 0.02  # Adjust this length as needed
    margin = 0.005

    # Calculate start and end points of the tick (perpendicular to the curve)
    x_tick_start = x_smooth[idx] - tick_length * perp_dx
    x_tick_end = x_smooth[idx] + tick_length * perp_dx
    y_tick_start = y_smooth[idx] - tick_length * perp_dy
    y_tick_end = y_smooth[idx] + tick_length * perp_dy

    if abs(perp_dx) < abs(perp_dy):

        ha = 'center'
        va = 'top'
        margin_x = 0
        margin_y = margin

    else:

        ha = 'right'
        va = 'center'

        margin_x = margin
        margin_y = 0

    axes[1].text(

        x_tick_end - margin_x, y_tick_end - margin_y,
        "{:0.2f}".format(omega), fontsize=12, ha=ha, va=va,
        color=c_smooth, # bbox=dict(facecolor='white', alpha=0.1)

    )

    # Plot the tick as a line
    axes[1].plot([x_tick_start, x_tick_end], [y_tick_start, y_tick_end], c=c_smooth, linewidth=2)

# axes[1].set_xlim((0.2, 0.8))

# <------------------ FINAL CALCULATIONS ------------------------------------------>
for ax in axes:
    ax.set_xlabel(x_label, fontsize=label_font_size, loc='right')
    ax.set_ylabel(y_label, fontsize=label_font_size, loc='top')
    ax.tick_params(axis='both', which='major', labelsize=label_font_size)

plt.tight_layout(pad=tight_pad)
plt.show()