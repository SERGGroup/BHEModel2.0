# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.constants import CALCULATION_FOLDER
from UNISIMConnect import UNISIMConnector
from tqdm import tqdm
import numpy as np
import os


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "4 - Water Based Case.usc")
output_folder = os.path.join(base_folder, "00 - Output", "4 - Overall Optimization Result - Water")

depth_cell = "C3"
grad_cell = "C4"
t_sep_perc_cell = "C8"

result_cells = {

    "M_RATIO": {"COMP": "F4", "HP": "F6"},
    "W_REL_NET": {"COMP": "F14", "HP": "F15"}

}

error_cell = "C14"


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #
n_geo = 10

# 1. CALCULATION RANGES -------------------------------------------------- >
depth_list = np.round(np.linspace(1000, 5000, n_geo), 1)
grad_list = np.round(np.linspace(20, 100, n_geo), 1)
t_sg_perc_list = np.linspace(0.1, 1, n_geo+1)[:-1]
depth_list, grad_list, t_sg_perc_list = np.meshgrid(depth_list, grad_list, t_sg_perc_list, indexing='ij')

# 3. RESULT ARRAY INIT -------------------------------------------------- >
new_shape = (depth_list.shape[0], depth_list.shape[1], depth_list.shape[2], 8)
output_array = np.empty(new_shape)
output_array[:] = np.nan

output_array[:, :, :, 0] = depth_list
output_array[:, :, :, 1] = grad_list
output_array[:, :, :, 2] = t_sg_perc_list

np.save(os.path.join(output_folder, 'output_array.npy'), output_array)

# 4. CALCULATIONS     -------------------------------------------------- >
pbar = tqdm(desc="Calculation Ongoing", total=n_geo**3)
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("OUTPUT RESULTS")

    # 4.1 CHANGE DEPTH ------------------------------------------------- >
    for i in range(new_shape[0]):

        spreadsheet.set_cell_value(depth_cell, output_array[i, 0, 0, 0])
        unisim.wait_solution()

        # 4.2 CHANGE GRADIENT  ---------------------------------------- >
        for j in range(new_shape[1]):

            spreadsheet.set_cell_value(grad_cell, output_array[i, j, 0, 1])
            unisim.wait_solution()

            # 4.3 CHANGE T_sep_%  ------------------------------------- >
            for k in range(new_shape[2]):

                spreadsheet.set_cell_value(t_sep_perc_cell, output_array[i, j, k, 2])
                unisim.wait_solution()

                # 4.4 COLLECT RESULTS  -------------------------------- >
                n = 3
                for system in ["COMP", "HP"]:

                    for velue in ["M_RATIO", "W_REL_NET"]:

                        output_array[i, j, k, n] = spreadsheet.get_cell_value(result_cells[velue][system])
                        n += 1

                output_array[i, j, k, 7] = spreadsheet.get_cell_value(error_cell)
                pbar.update(1)

                # 4.5 SAVE OUTPUT ARRAY  ------------------------------ >
                np.save(os.path.join(output_folder, 'output_array.npy'), output_array)

pbar.close()
