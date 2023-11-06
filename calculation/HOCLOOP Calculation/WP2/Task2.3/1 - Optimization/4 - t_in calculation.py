# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLRocksInfo,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Input data from:
#       "Optimization-T2.2&2.3.pptx" (in "0 - Resources" Folder)
#

# Fixed Values
p_in = 1.0          # [MPa]
t_surf = 15         # [C]
c_rock = 0.900      # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

tub_id = 0.0850     # [m]
tub_od = 0.1400     # [m]

# Variable Values (fixed for me)
mass_flow = 5       # [kg/s]
depth = 4000        # [m]
l_horiz = 5000      # [m]

# horizontal section specific
k_rock_horiz = 3.0       # [W/(m K)]
cas_od_horiz = 0.1940    # [m]
cas_id_horiz = 0.1750    # [m]

#vertical section specific
k_rock_vert = 2.0        # [W/(m K)]
cas_od_vert = 0.2730     # [m]
cas_id_vert = 0.2530     # [m]

# Variable Values (variable for me)
t_in_list = [25, 30, 35]    # [C]
k_ins = 0.0264              # [W/(m K)]
t_grad = 0.03               # [C/m]


# %%------------   CALCULATION                            -----------------------------------------------------------> #
for i in range(len(t_in_list)):

    t_in = t_in_list[i]
    t_rock = t_surf + t_grad * depth

    hs_geometry_vert = REELWELLGeometry(

        depth,
        tub_id=tub_id,
        tub_od=tub_od,
        cas_id=cas_id_vert,
        cas_od=cas_od_vert,
        k_insulation=k_ins,
        hot_in_tubing=True,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=False,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False,

        rocks_info=REELWELLRocksInfo(

            k_rocks=k_rock_vert

        )

    )

    hs_geometry_horiz = REELWELLGeometry(

        l_horiz,
        tub_id=tub_id,
        tub_od=tub_od,
        cas_id=cas_id_vert,
        cas_od=cas_od_vert,
        k_insulation=k_ins,
        hot_in_tubing=True,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=False,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False,

        rocks_info=REELWELLRocksInfo(

            k_rocks=k_rock_horiz

        )

    )


    bhe_in = PlantThermoPoint(["Water"], [1])
    bhe_in.set_variable("T", t_in)
    bhe_in.set_variable("P", p_in)
    bhe_in.m_dot = mass_flow

    well = REELWEELBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
        k_rocks=k_rock_vert, c_rocks=c_rock, rho_rocks=rho_rock,
        rw_geometry=hs_geometry_vert, max_iteration=20

    )

    heating_section = REELWELLHeatingSection(

        well, hs_geometry_horiz,
        integrate_temperature=True,

    )
    well.heating_section = heating_section

    main_time_points = [1, 7, 15, 30, 60, 90, 180, 365, 730, 1825, 3650, 7300]

    time_points = []
    time_points.extend(main_time_points)
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

    pbar = tqdm(desc="Calculate t_in={}".format(t_in), total=len(time_points))
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

    RES_FOLDER = os.path.join(

        constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
        "WP2", "Task2.3", "0 - Resources", "output", "optimization"

    )

    file_path = os.path.join(RES_FOLDER, "optimization t_in={:.0f}.xlsx".format(t_in))

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
