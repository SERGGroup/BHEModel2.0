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


BASE_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP3", "Task3.2", "EGC Paper"

)


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

#
#   Validation proposed by Pietro Ungar
#   Check Slides: ""
#
same_well_as_water = True

t_in = 10           # [C]
mass_flow = 30      # [kg/s]

if same_well_as_water:
    depth = 4250        # [m]
    l_overall = 8000    # [m]

else:
    depth = 3250        # [m]
    l_overall = 6000    # [m]

t_surf = 10         # [C]
t_grad = 0.05       # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.930      # [kJ/(kg K)]
rho_rock = 2360     # [kg/m^3]

t_rock = t_surf + t_grad * depth
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(

    l_horiz,
    tub_id=0.097,
    tub_od=0.143,
    cas_id=0.279,
    cas_od=0.297,
    hot_in_tubing=True,
    max_back_time=4,
    alpha_old=0.5,
    k_insulation=0.01,
    neglect_internal_heat_transfer=False,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in + 1)
bhe_in.set_variable("Q", 0.)
bhe_in.set_variable("P", bhe_in.get_variable("P"))
bhe_in.set_variable("T", t_in)
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

time_points = [7, 15, 30, 91, 182.5, 365, 730, 1825, 3650]

step = 50
profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]
profile_positions_vert = np.linspace(0, depth, int(depth / step + 1))[1:]
overall_profile_positions = np.concatenate((profile_positions_vert, profile_positions + depth))
x_pos = np.array(overall_profile_positions)
x_pos = np.reshape(x_pos, (x_pos.size, 1))

time_list = list()
t_out_list = list()
p_out_list = list()
w_out_list = list()
t_profile_list = list()
p_profile_list = list()
rho_profile_list = list()
h_profile_list = list()
rho_der_list = list()
h_der_list = list()

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

    rho_prof = np.concatenate((rho_list_vert.T, rho_list.T)).T
    h_prof = np.concatenate((h_list_vert.T, h_list.T)).T

    rho_profile_list.append(rho_prof)
    h_profile_list.append(h_prof)
    der_res = list()
    for prof in [rho_prof, h_prof]:

        prof_dx = (prof[:, 1:] - prof[:, :-1]) / (x_pos[1:] - x_pos[:-1]).T
        prof_dx[1, :] = -prof_dx[1, :]
        der_res.append(prof_dx)

    rho_der_list.append(der_res[0])
    h_der_list.append(der_res[1])

    pbar.update(1)

pbar.close()


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
other_name = ""

if hs_geometry.neglect_internal_heat_transfer:
    other_name += " - heat_negl"

if same_well_as_water:
    other_name += " - same_well"

file_path = os.path.join(BASE_FOLDER, "0 - Output", "co2_results - {}{}.xlsx".format(mass_flow, other_name))
data_exporter = {

    "well": well,
    "time_list": time_list,
    "t_out_list": t_out_list,
    "w_out_list": w_out_list,
    "p_out_list": p_out_list,
    "t_profile_list": t_profile_list,
    "p_profile_list": p_profile_list,
    "rho_derivative_list": rho_der_list,
    "h_derivative_list": h_der_list,
    "profile_positions": overall_profile_positions,

}

export_profiles_to_excel(file_path, data_exporter, reverse_time_position=False)
