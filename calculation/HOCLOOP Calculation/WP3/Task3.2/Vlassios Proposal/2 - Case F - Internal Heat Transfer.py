# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
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


BASE_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP3", "Task3.2", "Vlassios Proposal"

)


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

#
#   Validation proposed by Pietro Ungar
#   Check Slides: ""
#

t_in = 10           # [C]
p_in = 5            # [MPa]
mass_flow = 5       # [kg/s]

depth = 4500        # [m]
l_overall = 6500    # [m]

t_surf = 20         # [C]
t_grad = 0.04       # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.930      # [kJ/(kg K)]
rho_rock = 2360     # [kg/m^3]

t_rock = t_surf + t_grad * depth
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(

    l_horiz,
    tub_id=0.097,
    tub_od=0.143,
    cas_id=0.279,
    cas_od=0.297,
    hot_in_tubing=True,
    max_back_time=4,
    alpha_old=0.5,
    k_insulation=0.01,
    neglect_internal_heat_transfer=False,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)
hs_geometry.ra = 0.000045


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", p_in)
bhe_in.m_dot = mass_flow

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    rw_geometry=hs_geometry, max_iteration=20

)

heating_section = REELWELLHeatingSection(

    well, hs_geometry,
    integrate_temperature=True,

)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5

main_time_points = [7, 15, 30, 91, 182.5, 365, 730, 1825, 3650]

time_points = []
time_points.extend(main_time_points)
time_points.extend([7.78, 16.35, 28.41, 85.25, 176.99, 367.23, 761.72, 1895.72])
time_points.extend([0.08, 0.168, 0.25, 0.5, 1, 1.25, 2.5, 5])
time_points.sort()

step = 50
profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]
profile_positions_vert = np.linspace(0, depth, int(depth / step + 1))[1:]
overall_profile_positions = np.concatenate((profile_positions_vert, profile_positions + depth))
x_pos = np.array(overall_profile_positions)
x_pos = np.reshape(x_pos, (x_pos.size, 1))

time_list = list()
t_out_list = list()
p_out_list = list()
w_out_list = list()
t_profile_list = list()
p_profile_list = list()
rho_profile_list = list()
h_profile_list = list()
rho_der_list = list()
h_der_list = list()

pbar = tqdm(desc="Calculation Ongoing", total=len(time_points))
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

    rho_prof = np.concatenate((rho_list_vert.T, rho_list.T)).T
    h_prof = np.concatenate((h_list_vert.T, h_list.T)).T

    rho_profile_list.append(rho_prof)
    h_profile_list.append(h_prof)
    der_res = list()
    for prof in [rho_prof, h_prof]:

        prof_dx = (prof[:, 1:] - prof[:, :-1]) / (x_pos[1:] - x_pos[:-1]).T
        prof_dx[1, :] = -prof_dx[1, :]
        der_res.append(prof_dx)

    rho_der_list.append(der_res[0])
    h_der_list.append(der_res[1])

    pbar.update(1)

pbar.close()


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
file_path = os.path.join(BASE_FOLDER, "0 - Output", "case-2.xlsx")
data_exporter = {

    "well": well,
    "time_list": time_list,
    "t_out_list": t_out_list,
    "w_out_list": w_out_list,
    "p_out_list": p_out_list,
    "t_profile_list": t_profile_list,
    "p_profile_list": p_profile_list,
    "rho_derivative_list": rho_der_list,
    "h_derivative_list": h_der_list,
    "profile_positions": overall_profile_positions,

}

export_profiles_to_excel(file_path, data_exporter, reverse_time_position=False)


# %%------------   PLOT TIME VARIABLES                    -----------------------------------------------------------> #
partner_result_folder = os.path.join(BASE_FOLDER, "1 - Comparison", "partner results")
file_path_pressure = os.path.join(partner_result_folder, "pressure.csv")
file_path_temperature = os.path.join(partner_result_folder, "temperatures.csv")

data_pressure = np.loadtxt(file_path_pressure, delimiter=',', skiprows=2)
data_temperature = np.loadtxt(file_path_temperature, delimiter=',', skiprows=2)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4.5))
time_array = np.array(time_list)

partners_labels = ["IFPEN", "IFE"]
ax_titles = ["Pressure", "Temperature"]
ax_labels = ["Pressure [bar]", "Temperature [°C]"]

unifi_results = [np.array(p_out_list), np.array(t_out_list)]
partner_results = [data_pressure, data_temperature]
multiplier_UNIFI = [10, 1]
multiplier_partner = [100, 1]

time_corr = 24
#time_corr = 10

for n, ax in enumerate(axs):

    ax.scatter(time_array, unifi_results[n] * multiplier_UNIFI[n], marker="x", label="UNIFI", color="tab:green")
    partner_result = partner_results[n]

    for i in range(2):

        ax.plot(partner_result[:, 2*i] / time_corr, partner_result[:, 2*i + 1] * multiplier_partner[n], label=partners_labels[i])

    # ax.set_xscale("log")
    ax.legend()
    ax.set_title(ax_titles[n])
    ax.set_ylabel(ax_labels[n])
    ax.set_xlabel("time [days]")

plt.tight_layout()
plt.savefig(os.path.join(BASE_FOLDER, "1 - Comparison", "comparison.png"), dpi=300)
plt.show()


# %%------------   PLOT PROFILES                          -----------------------------------------------------------> #
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
colors = ["tab:blue", "tab:orange", "tab:green"]

for i, n in enumerate([0, 10, -1]):

    axs[0].plot(x_pos[1:], rho_der_list[n].T, colors[i])
    axs[1].plot(x_pos[1:], h_der_list[n].T, colors[i])

plt.tight_layout()
plt.show()