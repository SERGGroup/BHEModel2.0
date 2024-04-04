# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from pandas import DataFrame, ExcelWriter
from main_code import constants
from datetime import datetime
from tqdm import tqdm
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

mass_flow_base = 9.50149 # [kg/s]

cas_id_base = 0.1617     # [m]
cas_id_modifier = 1
cas_id = cas_id_base * cas_id_modifier

depth = 3000        # [m]
l_overall = 6500    # [m]

t_surf = 10         # [C]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_in = 10           # [C]
p_in = 3            # [MPa]
eta_HP = 0.4
eta_pump = 0.8

Le = 20
i_rate = 0.04
om_ratio = 0.05
hy = 8000
c_el = 0.075
alpha = (1 - (1 + i_rate) ** (-Le)) / i_rate
beta = (1 + alpha * om_ratio) / (alpha * hy)

bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)

s_0 = bhe_in.get_variable("s")
h_0 = bhe_in.get_variable("h")
T_0 = t_in + 273.15


# %%------------   INIT CALCULATIONS                      -----------------------------------------------------------> #
time = 3650

n_T_max = 10
n_grad_T = 6
n_m_dot = 9
n_t_in_water = 16

t_max_arr = np.linspace(80, 120, n_T_max)
grad_T_arr = np.linspace(50, 75, n_grad_T)
mass_flow_arr = np.linspace(4, 12, n_m_dot)
t_in_water = np.linspace(20, 80, n_t_in_water)

t_max_arr, grad_T_arr, mass_flow_arr, t_in_water = np.meshgrid(t_max_arr, grad_T_arr, mass_flow_arr, t_in_water, indexing='ij')
calculation_dict = {}
inputs_names = ["t_max", "grad_T", "m_dot_well", "T_in_well"]

for i in range(n_T_max):

    for j in range(n_grad_T):

        for k in range(n_m_dot):

            for n in range(n_t_in_water):

                key = "{}-{}-{}-{}".format(i, j, k, n)

                calculation_dict.update({

                    key:  [

                        t_max_arr[i, j, k, n],
                        grad_T_arr[i, j, k, n],
                        mass_flow_arr[i, j, k, n],
                        t_in_water[i, j, k, n]

                    ]
                })


# %%------------   CALCULATE                              -----------------------------------------------------------> #
pbar = tqdm(desc="Calculation Ongoing", total=len(calculation_dict.keys()))
results = {}

for key in calculation_dict.keys():

    t_max = calculation_dict[key][0]
    t_in = calculation_dict[key][3]

    bhe_in.set_variable("T", t_in)
    bhe_in.set_variable("P", p_in)

    mass_flow = calculation_dict[key][2]
    t_rock = t_surf + calculation_dict[key][1] / 1000 * depth
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
    results_list = []
    try:

        bhe_in.m_dot = mass_flow
        well.heating_section.time = time / 365
        well.update()

        if len(well.points) > 2 and well.points[-1].get_variable("rho") > 0:

            p_out = well.points[-1].get_variable("P")
            rho = well.points[-1].get_variable("rho")
            c_well = 1.15 * 1.05 * 2.86 * (0.105 * l_overall ** 2 + 1776 * l_overall * cas_id + 2.735E5)

            w_dot_hp = (1 - (t_in + 273.15)/(t_max + 273.15)) / eta_HP * well.power
            w_dot_hp = max(w_dot_hp, 0)
            c_HP = 0.33667 * w_dot_hp / 1000

            w_dot_pump = (p_in - p_out) * mass_flow / (rho * eta_pump) / 1000
            w_dot_pump = max(w_dot_pump, 0)

            w_dot_tot = (w_dot_hp + w_dot_pump)
            q_dh = w_dot_hp + well.power

            LCOH = ((c_HP + c_well) * beta + w_dot_tot * c_el) / q_dh
            results_list = [LCOH, q_dh, w_dot_tot, t_max]

    except:

        pass

    results.update({key: results_list})
    pbar.update(1)

pbar.close()
# %%
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "2 - EES calculation",

)
now = datetime.now()
now_str = "{}-{:02}-{:02} - {:02}.{:02}".format(now.year, now.month, now.day, now.hour, now.minute)
file_path = os.path.join(RES_FOLDER, "00 - Results", "calculation results files", "water-{}.xlsx".format(now_str))
df_results = list()

results_names = ["LCOH", "Q_DH", "w_net", "t_max_reach"]

for i in range(n_T_max):

    for j in range(n_grad_T):

        for k in range(n_m_dot):

            for n in range(n_t_in_water):

                key = "{}-{}-{}-{}".format(i, j, k, n)
                sub_list = calculation_dict[key]
                sub_list.extend(results[key])
                df_results.append(sub_list)

names = inputs_names
names.extend(results_names)
df = DataFrame(df_results)
df.columns = names
df.to_excel(file_path, index=False, startrow=1, startcol=1)
