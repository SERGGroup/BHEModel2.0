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
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#       "same diameter throughout the well" considered
#

p_in = 10           # [MPa]
t_in = 35           # [C]

mass_flow = 10      # [kg/s]

depth = 4212        # [m]
l_overall = 5500    # [m]

cas_id = 0.313614   # [m] (12.347 inch)
dd_cas = 0.015      # [m]
dd_ann = 0.0868     # [m]
dd_tub = 0.03       # [m]

t_grad = 0.020      # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_in + t_grad * depth - 10
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(

    l_horiz,
    tub_id=cas_id - dd_ann - dd_tub,
    tub_od=cas_id - dd_ann,
    cas_id=cas_id,
    cas_od=cas_id + dd_cas,
    k_insulation=0.01,
    hot_in_tubing=True,
    max_back_time=6,
    alpha_old=0.85,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_in,
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
main_time_points = [1, 2, 4, 8, 20]
time_points = np.array(main_time_points) * 365

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

bhe_in.m_dot = mass_flow
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
    "WP4", "5 - Comparison with Polish Calculation", "0 - Resources", "results"

)

file_path = os.path.join(RES_FOLDER, "result-p_in={}bar.xlsx".format(p_in*10))
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

ax.plot(time_array, p_out_list)
plt.xscale("log")
plt.show()
