# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from EESConnect import EESConnector
from main_code import constants
from datetime import datetime
from pandas import DataFrame
from tqdm import tqdm
import numpy as np
import os.path


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Geological data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

# Surface Plant Thermodynamic Parameters
dt_he = 40          # [C]
t_turb_in = 20      # [C]
eta_comp = 0.8      # [-]
eta_turb = 0.75     # [-]
t_max = 90          # [C]
t_sat = 10          # [C]

# Geological Data
grad_T = 50         # [C/km]
t_surf = 10         # [C]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

# Well design parameters
cas_id = 0.1617     # [m]
time = 3650         # [days]

# Optimization Parameters
l_hor_well_arr = np.linspace(500, 3500, 6)
depth_well_arr = np.linspace(500, 3500, 6)
m_dot_well_arr = np.linspace(5, 15, 11)

l_hor_well_ms, depth_well_ms, m_dot_well_ms = np.meshgrid(

    l_hor_well_arr, depth_well_arr, m_dot_well_arr, indexing='ij'

)

# Well input condition definition
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_sat)
bhe_in.set_variable("Q", 0)
p_sat = bhe_in.get_variable("P")

bhe_in.set_variable("T", t_sat)
bhe_in.set_variable("P", p_sat * 1.001)


# %%------------   EVALUATE WELL                          -----------------------------------------------------------> #
calculation_dict = {}
pbar = tqdm(desc="Well Calculation Ongoing", total=len(l_hor_well_arr)*len(depth_well_arr)*len(m_dot_well_arr))

for i in range(len(l_hor_well_arr)):

    for j in range(len(depth_well_arr)):

        for k in range(len(m_dot_well_arr)):

            l_horiz = l_hor_well_ms[i, j, k]
            depth = depth_well_ms[i, j, k]
            mass_flow = m_dot_well_ms[i, j, k]

            t_rock = t_surf + grad_T / 1000 * depth
            l_tot = l_horiz + depth

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

                    beta_well = well.points[-1].get_variable("P") / well.points[0].get_variable("P")
                    dt_well = well.points[-1].get_variable("T") - well.points[0].get_variable("T")

                    calculation_dict.update({

                        "{}-{}-{}".format(i, j, k): [

                            t_max, dt_he, t_turb_in,
                            eta_comp, eta_turb, dt_well,
                            beta_well, mass_flow, l_tot, t_rock

                        ]
                    })

            except:

                pass


            pbar.update(1)

pbar.close()


# %%------------   CALCULATE                              -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "2 - EES calculation",

)

EES_FILE = os.path.join(RES_FOLDER, "0 - EES Files", "base heat pump V2 - python.EES")

with EESConnector(EES_FILE, ees_decimal_separator=",", display_progress_bar=True) as ees:

    try:
        result = ees.calculate(calculation_dict)

    except:
        result = {}

print("Done!")

now = datetime.now()
now_str = "{}-{:02}-{:02} - {:02}.{:02}".format(now.year, now.month, now.day, now.hour, now.minute)
file_path = os.path.join(RES_FOLDER, "00 - Results", "calculation results files", "co2-opt-{}.xlsx".format(now_str))
df_results = list()

results_names = ["LCOH", "Q_DH", "w_net", "T_max", "P_max", "T_well", "P_well"]
inputs_names = [

    "t_max", "dt_he", "t_turb_in",
    "eta_comp", "eta_turb", "dT_well",
    "beta_well", "m_dot_well", "l_tot_well", "T_rocks"

]

for i in range(len(l_hor_well_arr)):

    for j in range(len(depth_well_arr)):

        for k in range(len(m_dot_well_arr)):

            key = "{}-{}-{}".format(i, j, k)
            if key in calculation_dict.keys():

                sub_list = calculation_dict[key]
                sub_list.extend(result[key])
                df_results.append(sub_list)

names = inputs_names
names.extend(results_names)
df = DataFrame(df_results)
df.columns = names
df.to_excel(file_path, index=False, startrow=1, startcol=1)
