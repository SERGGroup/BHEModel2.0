# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

mass_flow_base = 8.80149 # [kg/s]
mass_flow_modifier = 1
mass_flow = mass_flow_base * mass_flow_modifier

cas_id_base = 0.1617     # [m]
cas_id_modifier = 1
cas_id = cas_id_base * cas_id_modifier

depth = 750         # [m]
l_overall = 6500    # [m]

t_surf = 10         # [C]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_in = t_surf       # [C]
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0.)

p_sat = bhe_in.get_variable("P")
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_sat * 1.1)

s_0 = bhe_in.get_variable("s")
h_0 = bhe_in.get_variable("h")
T_0 = t_in + 273.15


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
time_points = [1/24, 1, 7, 15, 30, 90, 180, 365, ]
gradients = [0.050, 0.025]
pressure_losses = False

RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "0 - results",
    "Exergy Calculation", "Well Profiles",
    "CO2-Based"

)

pbar = tqdm(desc="Calculation Ongoing", total=len(time_points)*len(gradients))

for t_grad in gradients:

    t_rock = t_surf + t_grad * depth
    l_horiz = l_overall - depth
    hs_geometry = REELWELLGeometry(

        l_horiz,
        tub_id=0.1,
        tub_od=0.13,
        cas_id=cas_id,
        cas_od=cas_id + 0.015,
        k_insulation=0.1,
        hot_in_tubing=False,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=True,
        ignore_tubing_pressure_losses=(not pressure_losses),
        ignore_annulus_pressure_losses=(not pressure_losses)

    )

    step = 50
    profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]
    profile_positions_vert = np.linspace(0, depth, int(depth / step + 1))[1:]
    overall_profile_positions = np.concatenate((profile_positions_vert, profile_positions + depth))

    time_list = list()
    t_out_list = list()
    p_out_list = list()
    w_out_list = list()
    t_profile_list = list()
    h_profile_list = list()
    s_profile_list = list()
    ex_profile_list = list()

    bhe_in.m_dot = mass_flow
    for time in time_points:

        well = REELWEELBHE(

            bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
            k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
            rw_geometry=hs_geometry, max_iteration=30

        )

        heating_section = REELWELLHeatingSection(

            well, hs_geometry,
            integrate_temperature=True,

        )
        well.heating_section = heating_section
        well.heating_section.time = time / 365
        well.update()

        try:
            results_vert = well.get_property_profile(["h", "s", "T"], profile_positions_vert)
            results_hs = heating_section.get_property_profile(["h", "s", "T"], profile_positions)

        except:
            pass

        else:

            time_list.append(time)
            t_out_list.append(well.points[-1].get_variable("T"))
            p_out_list.append(well.points[-1].get_variable("P"))
            w_out_list.append(well.power)


            t_list_vert = results_vert[:, 2, :]
            h_list_vert = results_vert[:, 0, :] - 9.81 * profile_positions_vert.T / 1000 - h_0
            s_list_vert = results_vert[:, 1, :] - s_0
            ex_list_vert = h_list_vert - T_0 * s_list_vert

            t_list = results_hs[:, 2, :]
            h_list = results_hs[:, 0, :] - 9.81 * depth / 1000 - h_0
            s_list = results_hs[:, 1, :] - s_0
            ex_list = h_list - T_0 * s_list

            t_profile_list.append(np.concatenate((t_list_vert.T, t_list.T)).T)
            h_profile_list.append(np.concatenate((h_list_vert.T, h_list.T)).T)
            s_profile_list.append(np.concatenate((s_list_vert.T, s_list.T)).T)
            ex_profile_list.append(np.concatenate((ex_list_vert.T, ex_list.T)).T)

        pbar.update(1)

    file_path = os.path.join(

        RES_FOLDER, "pressure losses = {}".format(pressure_losses),
        "horizontal well - grad={}-{}m.xlsx".format(t_grad, depth)

    )

    data_exporter = {

        "well": well,
        "time_list": time_list,
        "t_out_list": t_out_list,
        "w_out_list": w_out_list,
        "p_out_list": p_out_list,
        "Temperature_profile_list": t_profile_list,
        "Enthalpy_profile_list": h_profile_list,
        "Entropy_profile_list": s_profile_list,
        "Exergy_profile_list": ex_profile_list,
        "profile_positions": overall_profile_positions

    }

    export_profiles_to_excel(file_path, data_exporter, reverse_time_position=False)

pbar.close()
