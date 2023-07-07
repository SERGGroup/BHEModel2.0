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
depth = 6000        # [m]
l_overall = 10000   # [m]
r_cur = 500         # [m]

t_rock = 172.5      # [C]
k_rock = 2.44       # [W/(m K)]
c_rock = 1.085      # [kJ/(kg K)]
rho_rock = 2542     # [kg/m^3]

t_surf = 15         # [°C]
geo_grad = 45       # [°C/km]

m_dot = 25          # [kg/s]

time_points = [0.1, 1, 10]


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(l_horiz)

bhe_in = PlantThermoPoint(["Carbon Dioxide", "Nitrogen"], [0.9, 0.1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("Q", 0)
tmp_point = bhe_in.duplicate()

well = SimplifiedBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    use_rk=True

)

heating_section = REELWELLHeatingSection(well, hs_geometry, hot_in_tubing=False, neglect_internal_heat_transfer=True)
well.heating_section = heating_section


# %%------------   OPTIMIZE DEPTH                         -----------------------------------------------------------> #
n_points = 20
depth_list = np.linspace(1500, 8000, n_points)

p_r_pump = 1
tmp_in_point = bhe_in.duplicate()
bhe_in.set_to_compression_result(tmp_in_point.get_variable("P") * p_r_pump, 0.8, tmp_in_point)

GEO_OPTIMUM_DICT = {}

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

    for depth in depth_list:

        well.dz_well = depth
        well.t_rocks = depth / 1e3 * geo_grad + t_surf
        hs_geometry.l_hor = l_overall - depth

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

        print("{} - {} -> {}, {}".format(depth, time, net_power, m_dot))

    GEO_OPTIMUM_DICT.update({

        time: {

            "depth": depth_list,
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

for time in time_points:

    ax_pow.plot(

        GEO_OPTIMUM_DICT[time]["depth"], GEO_OPTIMUM_DICT[time]["power"]["net"],
        label="{} years".format(time)

    )

ax_pow.set_xlabel(r'Well depth [m]', fontsize='large', loc='center')
ax_pow.set_ylabel(r'Power [kW]', fontsize='large', loc='center')
plt.legend()
plt.show()

# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

RES_FOLDER = os.path.join(constants.CALCULATION_FOLDER, "Initial HOCLOOP Analysis", "0 - Resources")
output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "6 - Depth Optimization.png"))
