# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.ndimage import gaussian_filter, uniform_filter
from main_code.constants import CALCULATION_FOLDER
from UNISIMConnect import UNISIMConnector
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import pandas as pd
import numpy as np
import random
import os


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "3 - Optimization-better.usc")
output_folder = os.path.join(base_folder, "00 - Output", "3 - Overall Optimization Results")

depth_cell = "C3"
grad_cell = "C4"
t_sg_perc_cell = "C12"

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

p_max_cell = "C14"
ihx_on_cell = "C15"
t_ihx_out_cell = "D15"
ihx_power_cell = "C23"


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


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #

# 1. CALCULATION RANGES -------------------------------------------------- >
depth_list = np.round(np.linspace(1000, 5000, 15), 1)
grad_list = np.round(np.linspace(30, 100, 15), 1)
X, Y = np.meshgrid(depth_list, grad_list, indexing='ij')

# 2. MINIMIZATION PARAM -------------------------------------------------- >
omega_values = np.linspace(0, 1, 21)[1:-1]
omega_list_calc = np.linspace(0, 1, 150)
omega_index = [np.argmin(np.abs(omega_list_calc - omega)) for omega in omega_values]

t_sg_min = 0.001
sep_perc_min = 0
n_calc = 30

# Refinement of results
n_fine = 150
method = 'cubic'

# Derivative Evaluation (Savgol Filter)
window_size_filter = 5
poly_order = 2

# Derivative Gaussian Smoothing
sigma_val = 0
der_limit = 4e-3

# Result Averaging
window_size = 5

t_sg_perc_list = np.linspace(t_sg_min, 1, n_calc)
t_sg_fine = np.linspace(np.min(t_sg_perc_list), np.max(t_sg_perc_list), n_fine)

sep_perc_list = np.linspace(sep_perc_min, 1 - sep_perc_min, n_fine)[:-1]
sep_perc_fine = sep_perc_list

t_sg_perc_list, sep_perc_list = np.meshgrid(t_sg_perc_list, sep_perc_list, indexing='ij')
t_sg_fine, sep_perc_fine = np.meshgrid(t_sg_fine, sep_perc_fine, indexing='ij')
calc_points = np.vstack((t_sg_perc_list.ravel(), sep_perc_list.ravel())).T

# 3. RESULT ARRAY INIT -------------------------------------------------- >
new_shape = (X.shape[0], X.shape[1], t_sg_perc_list.shape[0], t_sg_perc_list.shape[1])
t_sg_perc_max = np.empty(new_shape)
sep_perc_max = np.empty(new_shape)
m_dot_max = np.empty(new_shape)
w_rel_max = np.empty(new_shape)

p_max_list = np.empty(new_shape)
t_ihx_list = np.empty(new_shape)
ihx_power_list = np.empty(new_shape)

t_sg_perc_max[:] = np.nan
sep_perc_max[:] = np.nan
m_dot_max[:] = np.nan
w_rel_max[:] = np.nan

p_max_list[:] = np.nan
t_ihx_list[:] = np.nan
ihx_power_list[:] = np.nan

np.save(os.path.join(output_folder, 'depth_list.npy'), depth_list)
np.save(os.path.join(output_folder, 'grad_list.npy'), grad_list)

# 4. CALCULATIONS     -------------------------------------------------- >
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("CALCULATION")

    # 4.1 CHANGE DEPTH ------------------------------------------------- >
    for i in range(X.shape[0]):

        spreadsheet.set_cell_value(depth_cell, X[i, 0])
        unisim.wait_solution()

        # 4.2 CHANGE GRADIENT  ---------------------------------------- >
        for j in range(X.shape[1]):

            spreadsheet.set_cell_value(grad_cell, Y[i, j])
            unisim.wait_solution()

            # 4.3 EVALUATE DIFFERENT OPTIMIZATION PARAM --------------- >
            m_ratio_sublist = np.empty(t_sg_perc_list.shape)
            w_rel_sublist = np.empty(t_sg_perc_list.shape)
            m_ratio_sublist[:] = np.nan
            w_rel_sublist[:] = np.nan

            n = j + i * X.shape[1] + 1
            n_tot = X.shape[0] * X.shape[1]
            n_perc = n / n_tot * 100
            label_pbar = "calculating {} of {} ({:0.0f}%)".format(n, n_tot, n_perc)

            pbar = tqdm(desc=label_pbar, total=t_sg_perc_list.shape[0])
            for n in range(t_sg_perc_list.shape[0]):

                spreadsheet.set_cell_value(t_sg_perc_cell, t_sg_perc_list[n, 0])
                unisim.reset_iterators()
                timed_out = unisim.wait_solution(timeout=10)

                if not timed_out:

                    t_ihx_per_bound = np.array([0.2, 1])

                    for m in range(30):

                        t_ihx_per = np.mean(t_ihx_per_bound)
                        spreadsheet.set_cell_value(t_set_cell, t_ihx_per)
                        unisim.wait_solution(timeout=10)

                        try:

                            error = spreadsheet.get_cell_value(pinch_error_cell)
                            if abs(error) < 0.001:
                                break

                        except:

                            t_ihx_per_bound[0] = t_ihx_per

                        else:

                            if error > 0:
                                t_ihx_per_bound[0] = t_ihx_per

                            else:
                                t_ihx_per_bound[1] = t_ihx_per

                    try:
                        w_rel_sublist[n, :], m_ratio_sublist[n, :] = evaluate_params(

                            sheet=spreadsheet, sep_perc=sep_perc_fine[n, :], use_rel_ratio=True

                        )

                    except:

                        pass

                pbar.update(1)

            pbar.close()

            # 4.4 REFINE OPTIMIZATION PARAM GRID --------------- >
            nan_mask = np.where(~np.isnan(w_rel_sublist.ravel()))
            w_rel_fine = griddata(calc_points[nan_mask], w_rel_sublist.ravel()[nan_mask], (t_sg_fine, sep_perc_fine), method=method)
            m_ratio_fine = griddata(calc_points[nan_mask], m_ratio_sublist.ravel()[nan_mask], (t_sg_fine, sep_perc_fine), method=method)

            # 4.5 PREPARE MINIMIZATION           --------------- >
            a = np.log(w_rel_fine + 1)
            b = np.log(1 / m_ratio_fine)

            x_min_list = np.empty(omega_list_calc.shape)
            y_min_list = np.empty(omega_list_calc.shape)
            min_der_list = np.empty(omega_list_calc.shape)

            # 4.6 EVALUATE MIN FOR DIFFERENT OMEGA ------------ >
            for k, omega in enumerate(omega_list_calc):

                dx_smooth = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=1, delta=1)
                dy_smooth = savgol_filter((1 - omega) * a - omega * b, window_size_filter, poly_order, deriv=1, axis=0, delta=1)
                dtot = dx_smooth * dx_smooth + dy_smooth * dy_smooth
                dtot = gaussian_filter(dtot, sigma=sigma_val)

                min_der = np.nanmin(dtot)
                i_min = np.where(dtot == min_der)

                try:

                    x_min_list[k] = np.nanmean(t_sg_fine[i_min])
                    y_min_list[k] = np.nanmean(sep_perc_fine[i_min])
                    min_der_list[k] = min_der

                except:

                    pass

            # 4.7 FILTER OPTIMIZATION RESULTS ------------ >
            # Delete not minima
            min_der_list = (min_der_list - np.nanmin(min_der_list)) / (np.nanmax(min_der_list) - np.nanmin(min_der_list))
            not_minima = min_der_list > 4e-3
            x_min_list[not_minima] = np.nan
            y_min_list[not_minima] = np.nan

            # 4.8 SMOOTH OPTIMIZATION RESULTS ------------ >
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

            # 4.9 RETRIEVE FINAL RESULTS ------------ >
            for k, (idx, omega) in enumerate(zip(omega_index, omega_values)):

                # Append optimal T_SG_%
                spreadsheet.set_cell_value(t_sg_perc_cell, x_smooth[idx])
                unisim.wait_solution()

                w_rel_max[i, j, k], m_dot_max[i, j, k] = evaluate_params(

                    sheet=spreadsheet, sep_perc=y_smooth[idx], use_rel_ratio=True

                )

                # Append results
                t_sg_perc_max[i, j, k] = x_smooth[idx]
                sep_perc_max[i, j, k] = y_smooth[idx]

                p_max_list[i, j, k] = spreadsheet.get_cell_value(p_max_cell)
                t_ihx_list[i, j, k] = spreadsheet.get_cell_value(t_ihx_out_cell)
                ihx_power_list[i, j, k] = spreadsheet.get_cell_value(ihx_power_cell)

                np.save(os.path.join(output_folder, 'w_rel_max.npy'), w_rel_max)
                np.save(os.path.join(output_folder, 'm_dot_max.npy'), m_dot_max)
                np.save(os.path.join(output_folder, 't_sg_perc_max.npy'), t_sg_perc_max)
                np.save(os.path.join(output_folder, 'sep_perc_max.npy'), sep_perc_max)
                np.save(os.path.join(output_folder, 'p_max_list.npy'), p_max_list)
                np.save(os.path.join(output_folder, 't_ihx_list.npy'), t_ihx_list)
                np.save(os.path.join(output_folder, 'ihx_power_list.npy'), ihx_power_list)
