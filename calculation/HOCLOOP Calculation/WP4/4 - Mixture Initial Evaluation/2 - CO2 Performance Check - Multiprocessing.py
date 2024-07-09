from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
from datetime import datetime
import concurrent.futures
from tqdm import tqdm
import numpy as np
import time
import os

RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "4 - Mixture Initial Evaluation",
    "res"

)

support_folder = os.path.join(RES_FOLDER, "support")
base_filename = "res_{i}_{j}.npy"

# Surface Plant Thermodynamic Parameters
dt_he = 40          # [C]
t_turb_in = 20      # [C]
eta_comp = 0.8      # [-]
eta_turb = 0.75     # [-]
t_max = 90          # [C]
t_sat = 10          # [C]
other_comp = "Dimethyl ether"

# Geological Data
grad_T = 75         # [C/km]
t_surf = 10         # [C]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

# Well design parameters
cas_id = 0.1617     # [m]
l_horiz = 3600      # [m]
depth = 3000        # [m]
l_tot = l_horiz + depth

# Economic Parameter
Le = 20
i_rate = 0.04
om_ratio = 0.05
hy = 8000
c_el = 0.075
alpha = (1 - (1 + i_rate) ** (-Le)) / i_rate
beta = (1 + alpha * om_ratio) / (alpha * hy)

time_arr = [3.6, 36, 360, 3600]

n_grad_T = 7
n_m_dot = 11
n_conc = 11
n_time = len(time_arr)

grad_T_arr = np.linspace(35, 75, n_grad_T)
m_dot_well_arr = np.linspace(5, 15, n_m_dot)
conc_arr = np.linspace(0, 0.1, n_conc)

grad_T_arr, m_dot_well_arr, conc_arr, time_arr = np.meshgrid(grad_T_arr, m_dot_well_arr, conc_arr, time_arr, indexing='ij')

calculation_dict = {}

for i in range(n_grad_T):

    for j in range(n_m_dot):

        for k in range(n_conc):

            for n in range(n_time):

                calculation_dict.update({

                    "{}-{}-{}-{}".format(i, j, k, n):  [

                        grad_T_arr[i, j, k, n], m_dot_well_arr[i, j, k, n],
                        conc_arr[i, j, k, n], time_arr[i, j, k, n]

                    ]

                })

well_geometry = REELWELLGeometry(

    depth,
    tub_id=0.1,
    tub_od=0.13,
    cas_id=cas_id,
    cas_od=cas_id + 0.015,
    k_insulation=0.1,
    hot_in_tubing=False,
    max_back_time=4,
    alpha_old=0.5,
    neglect_internal_heat_transfer=True,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)

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

RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "4 - Mixture Initial Evaluation", "res",

)
pbar = tqdm(desc="Well Initialization Ongoing", total=len(calculation_dict.keys()))

for key in calculation_dict.keys():

    if len(calculation_dict[key]) <= 4:

        conc = calculation_dict[key][2]

        try:
            # Well input condition definition
            bhe_in = PlantThermoPoint(["Carbon Dioxide", other_comp], [1-conc, conc])
            bhe_in.set_variable("T", t_sat)
            bhe_in.set_variable("Q", 0)
            p_sat = bhe_in.get_variable("P")

            bhe_in.set_variable("T", t_sat)
            bhe_in.set_variable("P", p_sat * 1.001)

            calculation_dict[key].append(p_sat * 1.001)

        except:
            pass

    pbar.update(1)

pbar.close()

def evaluate_well(i, j):

    results = list()
    max_len = 0
    print("i={i} - j={j} - Calculation Started".format(i=i, j=j))

    for k in range(n_conc):

        for n in range(n_time):

            key = "{}-{}-{}-{}".format(i, j, k, n)
            grad_T = calculation_dict[key][0]
            mass_flow = calculation_dict[key][1]
            conc = calculation_dict[key][2]
            time = calculation_dict[key][3]

            sub_results = [grad_T, mass_flow, conc, time]

            t_rock = t_surf + grad_T / 1000 * depth

            try:

                # Well input condition definition
                bhe_in = PlantThermoPoint(["Carbon Dioxide", other_comp], [1 - conc, conc])
                bhe_in.set_variable("T", t_sat)
                bhe_in.set_variable("P", calculation_dict[key][4])
                bhe_in.m_dot = mass_flow

                well = REELWEELBHE(

                    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
                    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
                    rw_geometry=well_geometry, max_iteration=30

                )

                heating_section = REELWELLHeatingSection(

                    well, hs_geometry,
                    integrate_temperature=False,

                )
                well.heating_section = heating_section
                well.heating_section.time = time / 365
                well.update()

                if len(well.points) > 2 and well.points[-1].get_variable("rho") > 0:

                    print(well.points[-1].m_dot)
                    beta_well = well.points[-1].get_variable("P") / well.points[0].get_variable("P")
                    dt_well = well.points[-1].get_variable("T") - well.points[0].get_variable("T")
                    c_well = 1.15 * 1.05 * 2.86 * (0.105 * l_tot ** 2 + 1776 * l_tot * cas_id + 2.735E5)
                    LCOH = (c_well * beta) / well.power
                    LCOex = (c_well * beta) / (well.dex * mass_flow)

                    if dt_well > 0:

                        sub_results.append(beta_well)
                        sub_results.append(dt_well)

                        sub_results.append(well.dh)
                        sub_results.append(well.dex)
                        sub_results.append(well.dex_bottom)
                        sub_results.append(well.power)
                        sub_results.append(well.dex * mass_flow)
                        sub_results.append(LCOH)
                        sub_results.append(LCOex)

            except:

                pass

            if len(sub_results) > max_len:
                max_len = len(sub_results)

            results.append(sub_results)

    results_np = np.empty([len(results), max_len])
    results_np[:] = np.nan

    for k in range(len(results)):

        for n in range(len(results[k])):

            results_np[k, n] = results[k][n]

    curr_filename = os.path.join(support_folder, base_filename.format(i=i, j=j))
    np.save(curr_filename, results_np)

    print("i={i} - j={j} - Calculation Completed".format(i=i, j=j))


if __name__ == '__main__':

    print("Calculation Started!!")

    with concurrent.futures.ProcessPoolExecutor() as executor:

        future_list = list()

        for i in range(n_grad_T):

            for j in range(n_m_dot):

                time.sleep(0.5)
                future_list.append(executor.submit(evaluate_well, i, j))

    print("Calculation Completed!!")
    print("Saving Results ...")
    time.sleep(2)
    ovr_results = list()
    ovr_length = 0
    max_width = 0

    for i in range(n_grad_T):

        for j in range(n_m_dot):

            filename = os.path.join(support_folder, base_filename.format(i=i, j=j))

            if os.path.exists(filename):

                sub_arr = np.load(filename)
                ovr_length += sub_arr.shape[0]

                if max_width < sub_arr.shape[1]:
                    max_width = sub_arr.shape[1]

                ovr_results.append(sub_arr)

                # os.remove(filename)

    results = np.empty([ovr_length, max_width])
    results[:] = np.nan
    i = 0

    for sub_arr in ovr_results:

        for k in range(sub_arr.shape[0]):

            results[i, :] = sub_arr[k, :]
            i += 1

    now = datetime.now()
    now_str = "{}-{:02}-{:02} - {:02}.{:02}".format(now.year, now.month, now.day, now.hour, now.minute)
    filename = os.path.join(RES_FOLDER, "ovr_results-{}.npy".format(now_str))
    np.save(filename, results)

    print("Calculation Completed!!!")
