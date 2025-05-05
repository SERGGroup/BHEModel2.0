# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

#
#   Validation proposed by Pietro Ungar
#   Check Slides: ""
#

t_in = 10          # [C]
p_in = 5.5         # [MPa]
mass_flow = 5      # [kg/s]

cas_id = 0.279     # [m]

depth = 4500        # [m]
l_overall = 6500    # [m]

t_surf = 20         # [C]
t_grad = 0.04       # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.930      # [kJ/(kg K)]
rho_rock = 2360     # [kg/m^3]

t_rock = t_surf + t_grad * depth
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(

    l_horiz,
    tub_id=0.09715,
    tub_od=0.14339,
    cas_id=cas_id,
    cas_od=cas_id + 0.016,
    hot_in_tubing=True,
    max_back_time=4,
    alpha_old=0.5,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)

p_sat = bhe_in.get_variable("P")
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_sat * 1.001)
bhe_in.m_dot = mass_flow

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    rw_geometry=hs_geometry, max_iteration=20

)

heating_section = REELWELLHeatingSection(

    well, hs_geometry,
    integrate_temperature=True,

)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5

main_time_points = [7, 15, 30, 91, 182.5, 365, 730, 1825, 3650]

time_points = []
time_points.extend(main_time_points)
time_points.extend([7.78, 16.35, 28.41, 85.25, 176.99, 367.23, 761.72, 1895.72, 3931.18])
time_points.extend([0.08, 0.168, 0.25, 0.5, 1, 1.25, 2.5, 5])
time_points.sort()

step = 50
profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]
profile_positions_vert = np.linspace(0, depth, int(depth / step + 1))[1:]
overall_profile_positions = np.concatenate((profile_positions_vert, profile_positions + depth))

time_list = list()
t_out_list = list()
p_out_list = list()
w_out_list = list()
t_profile_list = list()
p_profile_list = list()

pbar = tqdm(desc="Calculation Ongoing", total=len(time_points))
for time in time_points:

    well.heating_section.time = time / 365
    well.update()

    time_list.append(time)
    t_out_list.append(well.points[-1].get_variable("T"))
    p_out_list.append(well.points[-1].get_variable("P"))
    w_out_list.append(well.power)

    t_list_vert, p_list_vert, rho_list_vert, h_list_vert = well.get_iteration_profile(profile_positions_vert)
    t_list, p_list, rho_list, h_list = heating_section.get_heating_section_profile(profile_positions)

    t_profile_list.append(np.concatenate((t_list_vert.T, t_list.T)).T)
    p_profile_list.append(np.concatenate((p_list_vert.T, p_list.T)).T)
    pbar.update(1)

pbar.close()


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP3", "Task3.2", "Vlassios Proposal", "0 - Output"

)

file_path = os.path.join(RES_FOLDER, "case-1.xlsx")
data_exporter = {

    "well": well,
    "time_list": time_list,
    "t_out_list": t_out_list,
    "w_out_list": w_out_list,
    "p_out_list": p_out_list,
    "t_profile_list": t_profile_list,
    "p_profile_list": p_profile_list,
    "profile_positions": overall_profile_positions

}

export_profiles_to_excel(file_path, data_exporter, reverse_time_position=False)


# %%------------   PLOT TIME VARIABLES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
time_array = np.array(time_list) / 365

ax.plot(time_array, np.array(p_out_list) * 10)
ax.plot(time_array, t_out_list)
plt.xscale("log")
# plt.ylim((70, 130))
plt.show()
