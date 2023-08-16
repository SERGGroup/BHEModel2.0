# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import REELWEELBHE, REELWELLGeometry
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-a.pdf" (in "0 - Resources\Case Studies" Folder)
#

t_in = 14.72        # [C]
depth = 1828.8      # [m]
mass_flow = 8.8     # [kg/s]

t_surf = 21.111     # [C]
t_grad = 0.01513    # [C/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_surf + t_grad * depth
hs_geometry = REELWELLGeometry(

    depth,
    tub_id=0.01,
    tub_od=0.011,
    cas_id=0.1617,
    cas_od=0.1778,
    hot_in_tubing=True,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=True,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 1)

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    t_surf=t_surf, rw_geometry=hs_geometry

)


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
profile_positions = np.linspace(0, depth, int(depth / step + 1))[1:]

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
    p_out_list.append(well.points[1].get_variable("P"))
    w_out_list.append(well.power)

    t_list, p_list, rho_list, h_list = well.get_iteration_profile(profile_positions)

    t_profile_list.append(t_list)
    p_profile_list.append(p_list)

    pbar.update(1)
    # print("{} -> {}".format(time, t_out_list[-1]))

pbar.close()


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP2", "Task2.2", "0 - Resources", "output"

)

file_path = os.path.join(RES_FOLDER, "case_a.xlsx")
data_exporter = {

    "well": well,
    "time_list": time_list,
    "t_out_list": t_out_list,
    "w_out_list": w_out_list,
    "p_out_list": p_out_list,
    "t_profile_list": t_profile_list,
    "p_profile_list": p_profile_list,
    "profile_positions": profile_positions

}

export_profiles_to_excel(file_path, data_exporter, reverse_time_position=False)


# %%------------   PLOT TIME VARIABLES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
time_array = np.array(time_list) / 365

ax.plot(time_array, t_out_list)
plt.xscale("log")
plt.show()
