# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from UNISIMConnect import UNISIMConnector
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np
import os

if os.name == "nt":
    from main_code.constants import CALCULATION_FOLDER

else:
    CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "old", "3 - general CO2 HTHP - downhole temperature + TSGperc.usc")
output_file = os.path.join(base_folder, "00 - Output", "2 - Test Different Temperature and Depth.csv")

depth_cell = "C3"
grad_cell = "C4"
sep_perc_cell = "C13"

m_ratio_cell = "G4"
m_ratio_rel_cell = "G6"

COP_cell = "C17"
eta_el_cell = "C19"
eta_cell = "C19"

T_sat_cell = "C22"

results = {

    depth_cell: list(),
    grad_cell: list(),
    sep_perc_cell: list(),

    m_ratio_cell: list(),
    m_ratio_rel_cell: list(),
    COP_cell: list(),
    eta_el_cell: list(),
    eta_cell: list(),

    T_sat_cell: list()

}

for i in range(9, 14, 1):

    results.update({"G{}".format(i): list()})

for i in range(16, 20, 1):

    results.update({"G{}".format(i): list()})


# %%------------   EVALUATE RESULTS                       -----------------------------------------------------------> #
depth_list = np.round(np.linspace(500, 5000, 15), 1)
grad_list = np.round(np.linspace(35, 100, 15), 1)

sep_perc_list = np.linspace(0, 1, 31)
sep_perc_list[0] = 0.001
sep_perc_list[-1] = 1 - 0.001

X, Y, Z = np.meshgrid(depth_list, grad_list, sep_perc_list, indexing='ij')

COP_list = np.empty(X.shape)
eta_list = np.empty(X.shape)
t_sat_list = np.empty(X.shape)
m_ratio_list = np.empty(X.shape)

COP_list[:] = np.nan
eta_list[:] = np.nan
t_sat_list[:] = np.nan
m_ratio_list[:] = np.nan

for key in results.keys():
    results[key] = list()

pbar = tqdm(desc="calculation", total=len(depth_list) * len(grad_list) * len(sep_perc_list))
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("CALCULATION")

    for i in range(X.shape[0]):

        spreadsheet.set_cell_value(depth_cell, X[i, 0, 0])
        unisim.wait_solution()

        for j in range(X.shape[1]):

            spreadsheet.set_cell_value(grad_cell, Y[i, j, 0])
            unisim.wait_solution()

            for k in range(X.shape[2]):

                spreadsheet.set_cell_value(sep_perc_cell, Z[i, j, k])
                unisim.wait_solution()

                COP_list[i, j, k] = spreadsheet.get_cell_value(COP_cell)
                eta_list[i, j, k] = spreadsheet.get_cell_value(eta_cell)
                t_sat_list[i, j, k] = spreadsheet.get_cell_value(T_sat_cell)
                m_ratio_list[i, j, k] = spreadsheet.get_cell_value(m_ratio_rel_cell)

                for key in results.keys():
                    results[key].append(spreadsheet.get_cell_value(key))

                pbar.update(1)

pbar.close()
pd.DataFrame(results).to_csv(output_file)


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
#specified_levels = np.linspace(0, 1, 50)
omega = 0.175
x_min = 1 / (20 * m_ratio_list) + omega * 3 / COP_list

plt.contourf(X[:, :, 0], Y[:, :, 0], (X*Y/1000 + 10)[:, :, 0], levels=20, cmap='viridis', extend='both')
plt.title('Contour Plot with Indexing')
plt.xlabel('Depth (m)')
plt.ylabel('Gradient (°C/m)')
plt.colorbar(label='Z value')
plt.show()