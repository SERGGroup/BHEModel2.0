# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.simplified_well.heating_sections.subclasses import REELWELLHeatingSection, REELWELLGeometry
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE
import matplotlib.pyplot as plt
from main_code import constants
import numpy as np
import os

# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

t_in = 30           # [C]
depth = 4500        # [m]
l_overall = 10000   # [m]
r_cur = 500         # [m]

t_rock = 172.5      # [C]
k_rock = 2.44       # [W/(m K)]
c_rock = 1.085      # [kJ/(kg K)]
rho_rock = 2542     # [kg/m^3]

m_dot = 25          # [kg/s]

l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(l_horiz)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #

bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0)

well = SimplifiedBHE(

    bhe_in, dz_well=depth, T_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    use_rk=True

)

heating_section = REELWELLHeatingSection(well, hs_geometry, hot_in_tubing=False, neglect_internal_heat_transfer=True)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #

n_points = 30
shape = 2.5
time_points = np.power((np.array(range(n_points)) + 1) / n_points, shape) * 10
RESULTS_DICT = {}

time_list = list()
t_out_list = list()
beta_out_list = list()
power_list = list()

bhe_in.m_dot = m_dot
tmp_point = bhe_in.duplicate()

for time in time_points:

    heating_section.time = time
    well.update()

    t_out = well.points[-1].get_variable("T")
    p_out = well.points[-1].get_variable("P")
    p_in = well.points[0].get_variable("P")

    tmp_point.set_to_expansion_result(p_in, 0.8, well.points[-1])
    power = (well.points[-1].get_variable("h") - tmp_point.get_variable("h")) * m_dot

    time_list.append(time)
    t_out_list.append(t_out)
    beta_out_list.append(p_in / p_out)
    power_list.append(power)

    print("{} - {} -> {}".format(m_dot, time, t_out_list[-1]))


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #

fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax = fig.add_subplot(1, 1, 1)
ax_pow = ax.twinx()

lines = list()
labels = list()

line_t = ax.plot(time_list, t_out_list, label="$T_{out}$")
line_p = ax_pow.plot(time_list, power_list, label="Power")

ax.set_xlabel(r'time [years]', fontsize='large', loc='center')
ax.set_ylabel(r'$T_{out}$ [°C]', fontsize='large', loc='center')
ax_pow.set_ylabel(r'Power [kW]', fontsize='large', loc='center')
plt.show()


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

RES_FOLDER = os.path.join(constants.CALCULATION_FOLDER, "Initial HOCLOOP Analysis", "0 - Resources")
output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "2 - Initial CO2 Evaluation.png"))
