# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import REELWEELBHE, REELWELLGeometry, REELWELLRocksInfo
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Data from:
#       "Deliverable 4.1 - Balmatt Case" (in "0 - Resources" Folder)
#

cas_id = 0.1617     # [m]

t_in = 65           # [C]
depth = 3600        # [m]
mass_flow = 8.8     # [kg/s]

t_surf = 10         # [C]
t_grad = 0.0325     # [C/m]
t_res_top = 130     # [C]
res_top = 3100      # [m]

k_rock = 2.68       # [W/(m K)]
c_rock = 0.93       # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_res_top + t_grad * (depth - res_top)

rocks_info = REELWELLRocksInfo(

    t_rocks=t_rock,
    k_rocks=k_rock,
    c_rocks=c_rock,
    rho_rocks=rho_rock,
    geo_gradient=t_grad,
    thermal_profile={

        0: 10,
        500: 25,
        900: 36,
        1600: 70,
        2200: 91,
        3100: 130,
        4500: 160

    }

)

hs_geometry = REELWELLGeometry(

    depth,
    tub_id=0.1,
    tub_od=0.13,
    cas_id=cas_id,
    cas_od=cas_id + 0.015,
    k_insulation=0.01,
    rocks_info=rocks_info,
    hot_in_tubing=True,
    max_back_time=3,
    alpha_old=0.5,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 1)

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    t_surf=t_surf, rw_geometry=hs_geometry, max_iteration=20

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5

main_time_points = [1, 7, 30, 90, 180, 365, 730, 1825, 3650, 7300]

time_points = []
time_points.extend(main_time_points)
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
    p_out_list.append(well.points[-1].get_variable("P"))
    w_out_list.append(well.power)

    t_list, p_list, rho_list, h_list = well.get_iteration_profile(profile_positions)

    t_profile_list.append(t_list)
    p_profile_list.append(p_list)

    pbar.update(1)
    #print("{} -> {}".format(time, t_out_list[-1]))
pbar.close()


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "1 - Balmatt Case", "0 - Resources", "results"

)

file_path = os.path.join(RES_FOLDER, "water-vertical.xlsx")
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

export_profiles_to_excel(file_path, data_exporter) #, times_in_main_tab=main_time_points)


# %%------------   PLOT TIME VARIABLES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
time_array = np.array(time_list) / 365

ax.plot(time_array, np.array(t_out_list))
ax.set_ylabel("T_out [°C]")
ax.set_xlabel("time [years]")
plt.xscale("log")
plt.show()
