# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import numpy as np


# %%------------   COMPARE WITH DISCRETIZATION            -----------------------------------------------------------> #
P_in = 5.74  # [MPa]
T_in = 20  # [°C]
depth = 7500  # [m]

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_in)
CO2_input.set_variable("P", P_in)

water_input = PlantThermoPoint(["Water"], [1])
water_input.set_variable("T", T_in)
water_input.set_variable("P", P_in)

bhe_CO2 = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=depth, t_rocks=400

)

bhe_CO2.update()
P_CO2_down = bhe_CO2.points[1].get_variable("P")
P_CO2_up = bhe_CO2.points[3].get_variable("P")

n_points = 1000
depth_list = np.linspace(0, depth, n_points, endpoint=True)
dz = depth_list[1] - depth_list[0]

down_point = bhe_CO2.points[0].duplicate()
up_point = bhe_CO2.points[3].duplicate()

p_up_list = list()
p_down_list = list()
rho_up_list = list()
rho_down_list = list()

for depth in depth_list[1:]:

    p_down = down_point.get_variable("P")
    h_down = down_point.get_variable("h")
    rho_down = down_point.get_variable("rho")

    p_up = up_point.get_variable("P")
    h_up = up_point.get_variable("h")
    rho_up = up_point.get_variable("rho")

    p_up_list.append(p_up)
    rho_up_list.append(rho_up)
    p_down_list.append(p_down)
    rho_down_list.append(rho_down)

    down_point.set_variable("P", p_down + rho_down * 9.81 * dz)
    down_point.set_variable("h", h_down + 9.81 * dz)

    up_point.set_variable("P", p_down - rho_down * 9.81 * dz)
    up_point.set_variable("h", h_down - 9.81 * dz)

plt.plot(p_down_list, rho_down_list)
plt.show()
