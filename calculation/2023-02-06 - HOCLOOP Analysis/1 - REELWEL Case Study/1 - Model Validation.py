# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.simplified_well.heating_sections.subclasses import REELWELLHeatingSection, REELWELLGeometry
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.other.error_band import draw_error_band
import matplotlib.pyplot as plt
from main_code import constants
import pandas as pd
import numpy as np
import os, math


# %%------------   IMPORT VALIDATION DATA                 -----------------------------------------------------------> #

#   Validation data from:
#       "20221107 HOCLOOP - CO2 case well" (in res folder)
#
#   Data extracted from plot using:
#       WebPlotDigitizer (https://automeris.io/WebPlotDigitizer)

RES_FOLDER = os.path.join(constants.CALCULATION_FOLDER, "2023-02-06 - HOCLOOP Analysis", "0 - Resources")
VALIDATION_FILE = os.path.join(RES_FOLDER, "validation_datasets", "t_over_time_validation_datasets.csv")
VALIDATION_ARR = np.array(pd.read_csv(VALIDATION_FILE).drop(0), dtype=float)
VALIDATION_DICT = {}
keys = [10, 5, 1]

for i in range(int(VALIDATION_ARR.shape[1] / 2)):

    j = 0

    for j in range(VALIDATION_ARR.shape[0]):

        if math.isnan(VALIDATION_ARR[j, 2 * i]):

            j -= 1
            break

    VALIDATION_DICT.update({

        keys[i]: {

            "t [years]": VALIDATION_ARR[:j+1, 2 * i],
            "T [°C]": VALIDATION_ARR[:j+1, 2 * i + 1]

        }

    })


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

t_in = 30           # [C]
depth = 4500        # [m]
l_overall = 10000   # [m]
r_cur = 500         # [m]

# !!! NOTE !!!
#
#   In the document provided by Ola (in the res folder) t_rock = 172°C while from the graphs used for the
#   validation (at pg. 2) such temperature appears to be around 160°C. This discrepancy, make me doubts about
#   the consistency of the other data.
#
#   For this reason, a correction factor has been introduced for the flow rate to overlap the calculated curve
#   and the validation points. The flow rate correction has been introduced only because it was easier to
#   implement, it simply reflects the fact that other parameters (such as rock conductivity or density) may
#   be different.
#
#   This has to be clarified with Ola.

t_rock = 160        # [C]
k_rock = 2.44       # [W/(m K)]
c_rock = 1.085      # [kJ/(kg K)]
rho_rock = 2542     # [kg/m^3]

l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(l_horiz)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #

bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 0.2)

well = SimplifiedBHE(

    bhe_in, dz_well=depth, T_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock

)

heating_section = REELWELLHeatingSection(well, hs_geometry, hot_in_tubing=False, neglect_internal_heat_transfer=True)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #

n_points = 30
shape = 2.5
time_points = np.power((np.array(range(n_points)) + 1) / n_points, shape) * 10
RESULTS_DICT = {}

for key in VALIDATION_DICT.keys():

    time_list = list()
    t_out_list = list()
    bhe_in.m_dot = key / 1.55

    for time in time_points:

        heating_section.time = time
        well.update()

        time_list.append(time)
        t_out_list.append(well.points[-1].get_variable("T"))

        print("{} - {} -> {}".format(key, time, t_out_list[-1]))

    RESULTS_DICT.update({

        key: {

            "t [years]": time_list,
            "T [°C]": t_out_list

        }

    })

# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #

fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax = fig.add_subplot(1, 1, 1)

lines = list()
labels = list()

for key in RESULTS_DICT.keys():

    x_calc = np.array(RESULTS_DICT[key]["t [years]"])
    y_calc = np.array(RESULTS_DICT[key]["T [°C]"])

    x_val = np.array(VALIDATION_DICT[key]["t [years]"])
    y_val = np.array(VALIDATION_DICT[key]["T [°C]"])

    line = draw_error_band(ax, x_calc, y_calc, 1.5, err_fixed=True, alpha=0.4, zorder=5)
    ax.scatter(x_val, y_val, s=80, marker="x", zorder=10, linewidths=1.5, color = line.get_color())

    lines.append(line)
    labels.append("$m_{{dot}}$ = {}[kg/s]".format(key))

ax.set_xlabel(r'time [years]', fontsize='large', loc='center')
ax.set_ylabel(r'$T_{out}$ [°C]', fontsize='large', loc='center')
plt.legend(lines, labels)
plt.show()


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "validation.png"))
