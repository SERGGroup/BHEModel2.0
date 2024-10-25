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
#

p_in = 1.4        # [MPa]
t_in = 10           # [C]

cas_id_base = 0.1617     # [m]
cas_id_modifier = 1
cas_id = cas_id_base * cas_id_modifier

t_surf = 11         # [C]
t_grad = 0.0325     # [C/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

eta_HP = 0.4
Le = 20
i_rate = 0.04
om_ratio = 0.05
hy = 8000
c_el = 0.075
alpha = (1 - (1 + i_rate) ** (-Le)) / i_rate
beta = (1 + alpha * om_ratio) / (alpha * hy)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0.)

p_in = bhe_in.get_variable("P") * 1.01
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
l_hor_well_arr = np.linspace(500, 3500, 6)
depth_well_arr = np.linspace(500, 3500, 6)
m_dot_well_arr = np.linspace(5, 15, 6)

w_out_list = list()
t_out_list = list()
p_out_list = list()
union = list()

pbar = tqdm(desc="Calculation Ongoing", total=(len(l_hor_well_arr))*(len(depth_well_arr))*(len(m_dot_well_arr)))
for i in range(len(m_dot_well_arr)):

    for j in range(len(depth_well_arr)):

        for k in range(len(l_hor_well_arr)):
            #print(m_dot_well_arr[i], l_hor_well_arr[j], depth_well_arr[k])
            union_list = list()
            l_horiz = l_hor_well_arr[k]
            union_list.append(l_horiz)

            depth = depth_well_arr[j]
            union_list.append(depth)

            mass_flow = m_dot_well_arr[i]
            union_list.append(mass_flow)

            union.append(union_list)
            l_overall=depth+l_horiz

            t_rock = t_surf + t_grad * depth
            step = 500
            profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]
            profile_positions_vert = np.linspace(0, depth, int(depth / step + 1))[1:]
            overall_profile_positions = np.concatenate((profile_positions_vert, profile_positions + depth))

            hs_geometry = REELWELLGeometry(

                l_horiz,
                tub_id=0.1,
                tub_od=0.13,
                cas_id=cas_id,
                cas_od=cas_id + 0.015,
                k_insulation=0.1,
                hot_in_tubing=True,
                max_back_time=4,
                alpha_old=0.5,
                neglect_internal_heat_transfer=False,
                ignore_tubing_pressure_losses=False,
                ignore_annulus_pressure_losses=False

            )

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

            well.update()
            bhe_in.m_dot = mass_flow
            t_out_list.append(well.points[-1].get_variable("T"))
            p_out_list.append(well.points[-1].get_variable("P"))
            t_list_vert, p_list_vert, rho_list_vert, h_list_vert = well.get_iteration_profile(profile_positions_vert)
            t_list, p_list, rho_list, h_list = heating_section.get_heating_section_profile(profile_positions)
            w_out_list.append(well.power)
            #print(well.power)
            print(t_list_vert)

            #calcolo del LCOH

            c_well = 1.15 * 1.05 * 2.86 * (0.105 * l_overall ** 2 + 1776 * l_overall * cas_id + 2.735E5)
            LCOH = (c_well * beta) / well.power

            print(LCOH)


            pbar.update(1)

print(union, w_out_list)