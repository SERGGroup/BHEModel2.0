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
mass_flow = 8.8       # [kg/s]

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
    k_insulation=0.1,
    rocks_info=rocks_info,
    hot_in_tubing=True,
    max_back_time=3,
    alpha_old=0.5,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "1 - Balmatt Case", "0 - Resources", "results"

)

shape = 2.5
n_points = 30
step = 50

rho_in_list = [250, 500, 750]
time_points = [1, 7, 30, 90, 180, 365, 730, 1825, 3650, 7300]
profile_positions = np.linspace(0, depth, int(depth / step + 1))[1:]
pbar = tqdm(desc="Calculation Ongoing", total=len(time_points)*len(rho_in_list))

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,6))

for rho_in in rho_in_list:

    bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
    bhe_in.set_variable("T", t_in)
    bhe_in.set_variable("rho", rho_in)

    well = REELWEELBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock,
        k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
        t_surf=t_surf, rw_geometry=hs_geometry, max_iteration=20

    )

    time_list = list()
    t_out_list = list()
    beta_list = list()
    w_out_list = list()
    t_profile_list = list()
    p_profile_list = list()

    bhe_in.m_dot = mass_flow
    for time in time_points:

        try:

            well.heating_section.time = time / 365
            well.update()

            time_list.append(time)
            t_out_list.append(well.points[-1].get_variable("T"))
            beta_list.append(well.points[-1].get_variable("P") / well.points[0].get_variable("P"))
            w_out_list.append(well.power)

            t_list, p_list, rho_list, h_list = well.get_iteration_profile(profile_positions)

            t_profile_list.append(t_list)
            p_profile_list.append(p_list)

        except:

            t_out_list.append(-1)
            beta_list.append(-1)
            w_out_list.append(-1)
            t_profile_list.append(list())
            p_profile_list.append(list())

        pbar.update(1)

    time_array = np.array(time_list) / 365
    axs[0].plot(time_array, np.array(t_out_list), label="rho_in = {}".format(rho_in))
    axs[1].plot(time_array, np.array(beta_list), label="rho_in = {}".format(rho_in))
    axs[2].plot(time_array, np.array(w_out_list), label="rho_in = {}".format(rho_in))

    file_path = os.path.join(RES_FOLDER, "CO2-vertical-{}.xlsx".format(rho_in))
    data_exporter = {

        "well": well,
        "time_list": time_list,
        "t_out_list": t_out_list,
        "w_out_list": w_out_list,
        "p_out_list": beta_list,
        "t_profile_list": t_profile_list,
        "p_profile_list": p_profile_list,
        "profile_positions": profile_positions

    }

    export_profiles_to_excel(file_path, data_exporter)  # , times_in_main_tab=main_time_points)

pbar.close()
plt.legend()
axs[0].set_ylabel("T_out [°C]")
axs[1].set_ylabel("Beta [-]")
axs[2].set_ylabel("W_out [kW]")

for ax in axs:

    ax.legend()
    ax.set_xscale("log")
    ax.set_xlabel("time [years]")


plt.tight_layout(pad=2)
plt.show()
