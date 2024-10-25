# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.constants import CALCULATION_FOLDER
from UNISIMConnect import UNISIMConnector
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np
import os


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
unisim_path = os.path.join(base_folder, "0 - UNISIM Files", "old", "3 - general CO2 HTHP - downhole temperature + TSGperc.usc")
output_file = os.path.join(base_folder, "00 - Output", "0.1 - first test.csv")

t_sg_perc_cell = "C14"
sep_perc_cell = "C13"

m_ratio_cell = "G4"
m_ratio_rel_cell = "G6"
COP_cell = "C17"
eta_el_cell = "C19"
eta_cell = "C19"

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
t_sg_perc_list = np.linspace(0, 1, 31)
t_sg_perc_list[0] = 0.001

sep_perc_list = np.linspace(0, 1, 31)
sep_perc_list[0] = 0.001
sep_perc_list[-1] = 1 - 0.001

X, Y = np.meshgrid(t_sg_perc_list, sep_perc_list, indexing='ij')

COP_list = np.empty(X.shape)
eta_list = np.empty(X.shape)
m_ratio_list = np.empty(X.shape)

COP_list[:] = np.nan
eta_list[:] = np.nan
m_ratio_list[:] = np.nan

for key in results.keys():
    results[key] = list()

pbar = tqdm(desc="calculation", total=len(t_sg_perc_list) * len(sep_perc_list))
with (UNISIMConnector(unisim_path, close_on_completion=False) as unisim):

    spreadsheet = unisim.get_spreadsheet("CALCULATION")

    for i in range(X.shape[0]):

        spreadsheet.set_cell_value(t_sg_perc_cell, X[i, 0])
        unisim.wait_solution()

        for j in range(X.shape[1]):

            spreadsheet.set_cell_value(sep_perc_cell, Y[i, j])
            unisim.wait_solution()
            COP_list[i, j] = spreadsheet.get_cell_value(COP_cell)
            m_ratio_list[i, j] = spreadsheet.get_cell_value(m_ratio_rel_cell)
            eta_list[i, j] = spreadsheet.get_cell_value(eta_cell)

            for key in results.keys():
                results[key].append(spreadsheet.get_cell_value(key))

            pbar.update(1)

pbar.close()
pd.DataFrame(results).to_csv(output_file)


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
#specified_levels = np.linspace(0, 1, 50)
omega = 0.175
x_min = 1 / (20 * m_ratio_list) + omega * 3 / COP_list

plt.contourf(X, Y, 1/x_min, levels=20, cmap='viridis', extend='both')
plt.title('Contour Plot with Indexing')
plt.xlabel('t_sg_perc')
plt.ylabel('sep_perc')
plt.colorbar(label='Z value')
plt.show()

