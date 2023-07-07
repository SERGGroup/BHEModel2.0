# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import REELWELLHeatingSection, REELWELLGeometry
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.optimize import minimize_scalar
from main_code import constants
import matplotlib.pyplot as plt
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
time_points = [0.1, 1, 10]
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(l_horiz, hot_in_tubing=False)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Carbon Dioxide"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0)
tmp_point = bhe_in.duplicate()

well = SimplifiedBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    use_rk=True

)

heating_section = REELWELLHeatingSection(well, hs_geometry, neglect_internal_heat_transfer=True)
well.heating_section = heating_section


# %%------------   OPTIMIZE INLET PRESSURE                -----------------------------------------------------------> #
n_points = 10
p_r_pump_list = np.linspace(1, 1.5, n_points)
tmp_in_point = bhe_in.duplicate()

P_IN_OPTIMUM_DICT = {}

for time in time_points:

    heating_section.time = time
    m_dot_list = list()
    power_dict = {

        "net": list(),
        "turb": list(),
        "pump": list(),

    }
    points_dict = {

        "bhe_in": list(),
        "hs_in": list(),
        "hs_out": list(),
        "bhe_out": list(),
        "turb_out": list(),
        "pump_in": list()

    }

    for p_r_pump in p_r_pump_list:

        bhe_in.set_to_compression_result(tmp_in_point.get_variable("P") * p_r_pump, 0.8, tmp_in_point)

        def opt_fun_m_dot(m_dot):

            net_power, turb_power, pump_power, tmp_point  = evaluate_system(m_dot)
            return -net_power

        def evaluate_system(m_dot):

            bhe_in.m_dot = m_dot
            pump_power = (bhe_in.get_variable("h") - tmp_in_point.get_variable("h")) * m_dot

            well.update()

            p_in = tmp_in_point.get_variable("P")
            tmp_point.set_to_expansion_result(p_in, 0.8, well.points[-1])
            turb_power = (well.points[-1].get_variable("h") - tmp_point.get_variable("h")) * m_dot
            net_power = turb_power - pump_power

            return net_power, turb_power, pump_power, tmp_point

        m_dot = minimize_scalar(opt_fun_m_dot, method='Bounded', bounds=(1, 50)).x
        net_power, turb_power, pump_power, tmp_point = evaluate_system(m_dot)

        m_dot_list.append(m_dot)
        power_dict["net"].append(net_power)
        power_dict["turb"].append(turb_power)
        power_dict["pump"].append(pump_power)

        points_dict["bhe_in"].append(well.points[0].duplicate())
        points_dict["hs_in"].append(well.points[1].duplicate())
        points_dict["hs_out"].append(well.points[2].duplicate())
        points_dict["bhe_out"].append(well.points[3].duplicate())
        points_dict["turb_out"].append(tmp_point.duplicate())
        points_dict["pump_in"].append(tmp_in_point.duplicate())

        print("{} - {} -> {}, {}".format(p_r_pump, time, net_power, m_dot))

    P_IN_OPTIMUM_DICT.update({

        time: {

            "pr_pump": p_r_pump_list,
            "m_dot": m_dot_list,
            "power": power_dict,
            "points": points_dict

        }

    })


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax_pow = fig.add_subplot(1, 1, 1)

time = time_points[0]

for key in P_IN_OPTIMUM_DICT[time]["power"].keys():

    ax_pow.plot(

        P_IN_OPTIMUM_DICT[time]["pr_pump"], P_IN_OPTIMUM_DICT[time]["power"][key],
        label="{} power".format(key)

    )

ax_pow.set_xlabel(r'$\beta_{pump}$ [-]', fontsize='large', loc='center')
ax_pow.set_ylabel(r'Power [kW]', fontsize='large', loc='center')
plt.legend()
plt.show()

# %%------------   PLOT POINTS                            -----------------------------------------------------------> #

fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax_pow = fig.add_subplot(1, 1, 1)

time = time_points[0]
var = "P"

for key in P_IN_OPTIMUM_DICT[time]["points"].keys():

    y_list = list()

    for point in P_IN_OPTIMUM_DICT[time]["points"][key]:

        y_list.append(point.get_variable(var))


    ax_pow.plot(

        P_IN_OPTIMUM_DICT[time]["pr_pump"], y_list,
        label="{}".format(key)

    )

ax_pow.set_xlabel(r'$\beta_{pump}$ [-]', fontsize='large', loc='center')
ax_pow.set_ylabel('{} [{}]'.format(var, bhe_in.get_unit(var)), fontsize='large', loc='center')
plt.legend()
plt.show()


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

RES_FOLDER = os.path.join(constants.CALCULATION_FOLDER, "Initial HOCLOOP Analysis", "0 - Resources")
output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "4 - Optimization CO2 pump power.png"))
