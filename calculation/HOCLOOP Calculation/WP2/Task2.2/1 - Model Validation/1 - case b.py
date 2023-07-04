# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.simplified_well.heating_sections.subclasses import REELWELLHeatingSection, REELWELLGeometry
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE
from main_code import constants
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

t_in = 45           # [C]
depth = 3000        # [m]
l_overall = 6500    # [m]
mass_flow = 8.80149 # [kg/s]

t_surf = 11         # [C]
t_grad = 0.0325     # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_surf + t_grad * depth
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(l_horiz)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 2.3)

well = SimplifiedBHE(

    bhe_in, dz_well=depth, T_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock

)

heating_section = REELWELLHeatingSection(well, hs_geometry, hot_in_tubing=True, neglect_internal_heat_transfer=True)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5

time_points = [1, 180, 365, 730, 1460, 2555, 3650]
time_points.extend([7.78, 16.35, 28.41, 85.25, 176.99, 367.23, 761.72, 1895.72, 3931.18])
time_points.extend([0.08, 0.168, 0.25, 0.5, 1, 1.25, 2.5, 5])
time_points.sort()

step = 500
profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]

time_list = list()
t_out_list = list()
p_out_list = list()
w_out_list = list()
profile_list = list()
bhe_in.m_dot = mass_flow

for time in time_points:

    heating_section.time = time / 365
    well.update()

    time_list.append(time)
    t_out_list.append(well.points[-1].get_variable("T"))
    p_out_list.append(well.points[-1].get_variable("P"))
    w_out_list.append(well.power)
    profile_list.append(heating_section.get_heating_section_profile(profile_positions))

    print("{} -> {}".format(time, t_out_list[-1]))


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP2", "Task2.2", "0 - Resources", "output"

)

file_path = os.path.join(RES_FOLDER, "case_b.xlsx")

main_data = {

    'Time': time_list,
    'T_out': t_out_list,
    'W_out': w_out_list

}
main_dataframe = pd.DataFrame(main_data)

n_profile_points = len(profile_list[0][0,:])
profile_list = np.array(profile_list)
profile_array = np.zeros((len(t_out_list), 2 * n_profile_points + 1))
profile_array[:, 0] = np.array(time_list)
profile_array[:, 1:(1+n_profile_points)] = profile_list[:, 0, :]
profile_array[:, (1+n_profile_points):] = profile_list[:, 1, :]

profile_dataframe = pd.DataFrame(profile_array)

writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
main_dataframe.to_excel(writer, sheet_name='Main Results', index=False)
profile_dataframe.to_excel(writer, sheet_name='Temperature Profiles', index=False)
writer.save()


# %%------------   PLOT PROFILES                          -----------------------------------------------------------> #
position_array = np.array(profile_positions)
position_array = np.concatenate((position_array, np.flip(position_array)))

for i in range(len(time_list)):

    plt.plot(profile_array[i, 1:], position_array)

plt.ylabel("horizontal distance [m]")
plt.xlabel("Temperature [°C]")
plt.show()


# %%------------   PLOT TIME VARIABLES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
time_array = np.array(time_list) / 365

ax.plot(time_array, t_out_list)
plt.xscale("log")
plt.show()
