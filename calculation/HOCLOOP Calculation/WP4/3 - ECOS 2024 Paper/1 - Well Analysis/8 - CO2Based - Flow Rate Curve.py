# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from pandas import DataFrame, ExcelWriter
from main_code import constants
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

mass_flow_base = 8.80149 # [kg/s]
cas_id_base = 0.1617     # [m]
cas_id_modifier = 1
cas_id = cas_id_base * cas_id_modifier

depth = 1000        # [m]
l_overall = 5500    # [m]

t_surf = 11         # [C]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_in = 10           # [C]
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0.)

p_sat = bhe_in.get_variable("P")
p_in = p_sat * 1.01
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)
s_0 = bhe_in.get_variable("s")
h_0 = bhe_in.get_variable("h")
T_0 = t_in + 273.15


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
mass_flow_modifiers = np.linspace(0.4, 1.5, 35)
gradients = [0.075, 0.05, 0.035]
time = 3650

results = np.zeros([len(gradients), len(mass_flow_modifiers), 10])
pbar = tqdm(desc="Calculation Ongoing", total=len(gradients)*len(mass_flow_modifiers))

for i in range(len(gradients)):

    for j in range(len(mass_flow_modifiers)):

        mass_flow_modifier = mass_flow_modifiers[j]
        mass_flow = mass_flow_base * mass_flow_modifier

        results[i, j, 0] = mass_flow
        for k in range(9):
            results[i, j, k + 1] = np.nan

        t_rock = t_surf + gradients[i] * depth
        l_horiz = l_overall - depth
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
            neglect_internal_heat_transfer=True,
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
            integrate_temperature=False,

        )
        well.heating_section = heating_section

        try:

            bhe_in.m_dot = mass_flow
            well.heating_section.time = time / 365
            well.update()

            if len(well.points) > 2 and well.points[-1].get_variable("rho") > 0:

                beta = well.points[-1].get_variable("P") / well.points[0].get_variable("P")
                dT = well.points[-1].get_variable("T") - well.points[0].get_variable("T")
                dh = well.points[-1].get_variable("h") - well.points[0].get_variable("h")
                ds = well.points[-1].get_variable("s") - well.points[0].get_variable("s")
                dex = dh - T_0 * ds

                results[i, j, 1] = well.power
                results[i, j, 2] = well.power * dex / dh
                results[i, j, 3] = dh
                results[i, j, 4] = dex
                results[i, j, 5] = beta
                results[i, j, 6] = dT
                results[i, j, 7] = well.points[-1].get_variable("T")
                results[i, j, 8] = well.points[-1].get_variable("P")
                results[i, j, 9] = well.points[-1].get_variable("rho")

        except:

            pass

        pbar.update(1)

pbar.close()


# %%------------  EXPORT DATA                             -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "0 - results",
    "Well Calculations", "Well Behaviour"

)

file_path = os.path.join(RES_FOLDER, "CO2-{}m-{}days.xlsx".format(depth, time))
col_names = ["flow rate", "power", "exergy", "dh", "dex", "beta", "dt", "T_out", "P_out", "rho_out"]

with ExcelWriter(file_path) as excel_writer:

    for i in range(len(gradients)):

        df = DataFrame(results[i, :, :])
        df.columns = col_names
        df.to_excel(excel_writer, sheet_name="grad={}".format(gradients[i]*1000))
