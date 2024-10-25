# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.constants import CALCULATION_FOLDER
from UNISIMConnect import UNISIMConnector
from tqdm import tqdm
import numpy as np
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
    "q_steam":      "H17",
    "ihx_power":    "C23",

}

m_ratio_cells = {

    "m_ratio":      "H4",
    "m_ratio_rel":  "H6",

}

p_max_cell = "C14"
ihx_on_cell = "C15"
t_ihx_out_cell = "D15"


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

    ihx_power = (1 - sep_perc) * sheet.get_cell_value(rel_power_cells["ihx_power"])

    return net_power, m_ratio, ihx_power


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #

# 1. CALCULATION RANGES -------------------------------------------------- >
depth_list = np.round(np.linspace(1000, 5000, 15), 1)
grad_list = np.round(np.linspace(30, 100, 15), 1)
X, Y = np.meshgrid(depth_list, grad_list, indexing='ij')

t_sg_min = 0.001
sep_perc_min = 0
n_calc = 30
n_fine = 150

t_sg_perc_list = np.linspace(t_sg_min, 1, n_calc)
sep_perc_list = np.linspace(sep_perc_min, 1 - sep_perc_min, n_fine)[:-1]
ones_array = np.ones(sep_perc_list.shape)

# 3. RESULT ARRAY INIT -------------------------------------------------- >
new_shape = (X.shape[0], X.shape[1], len(t_sg_perc_list), len(sep_perc_list))

m_ratio_list = np.empty(new_shape)
w_rel_list = np.empty(new_shape)
p_max_list = np.empty(new_shape)
t_ihx_list = np.empty(new_shape)
ihx_power_list = np.empty(new_shape)

m_ratio_list[:] = np.nan
w_rel_list[:] = np.nan
p_max_list[:] = np.nan
t_ihx_list[:] = np.nan
ihx_power_list[:] = np.nan

np.save(os.path.join(output_folder, 'depth_list.npy'), depth_list)
np.save(os.path.join(output_folder, 'grad_list.npy'), grad_list)
np.save(os.path.join(output_folder, 't_sg_perc_list.npy'), t_sg_perc_list)
np.save(os.path.join(output_folder, 'sep_perc_list.npy'), sep_perc_list)

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
            n = j + i * X.shape[1] + 1
            n_tot = X.shape[0] * X.shape[1]
            n_perc = n / n_tot * 100
            label_pbar = "calculating {} of {} ({:0.0f}%)".format(n, n_tot, n_perc)

            pbar = tqdm(desc=label_pbar, total=t_sg_perc_list.shape[0])
            for n in range(len(t_sg_perc_list)):

                spreadsheet.set_cell_value(t_sg_perc_cell, t_sg_perc_list[n])
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

                        w_rel_list[i, j, n, :], m_ratio_list[i, j, n, :], t_ihx_list[i, j, n, :] = evaluate_params(

                            sheet=spreadsheet, sep_perc=sep_perc_list, use_rel_ratio=True

                        )

                        p_max_list[i, j, n, :] = spreadsheet.get_cell_value(p_max_cell) * ones_array
                        t_ihx_list[i, j, n, :] = spreadsheet.get_cell_value(t_ihx_out_cell) * ones_array

                    except:

                        pass

                pbar.update(1)

            np.save(os.path.join(output_folder, 'w_rel_max.npy'), w_rel_list)
            np.save(os.path.join(output_folder, 'm_dot_max.npy'), m_ratio_list)
            np.save(os.path.join(output_folder, 'p_max_list.npy'), p_max_list)
            np.save(os.path.join(output_folder, 't_ihx_list.npy'), t_ihx_list)
            np.save(os.path.join(output_folder, 'ihx_power_list.npy'), ihx_power_list)

            pbar.close()
