# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.power_plants.HTHP.subclasses.CO2_heat_pump_py import CO2HeatPumpThermo
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.optimize import minimize
import numpy as np
import gc


# %%-------------------------------------   INITIALIZE CALCULATIONS             -------------------------------------> #
t_in_well = 50
grad = 50
depth = 1000
p_steam = 1
omega = 1

fluid = "Carbon Dioxide"
bhe_in = PlantThermoPoint([fluid], [1])
bhe_in.set_variable("T", bhe_in.RPHandler.TC)
bhe_in.set_variable("P", bhe_in.RPHandler.PC)
rho_crit = bhe_in.get_variable("rho")

bhe_in.set_variable("T", t_in_well)
bhe_in.set_variable("rho", rho_crit)
p_ref = bhe_in.get_variable("P")
q_ref = bhe_in.get_variable("Q")

t_rocks = t_in_well + grad * depth / 1000


def update_values(depth_in, t_rock, t_amb, p_in, T_SG_perc_curr):

    new_hthp = CO2HeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth_in, T_rock=t_rock,
        T_in_BHE=t_amb, P_in_BHE=p_in,
        T_ambient=t_amb - 10

    )

    new_hthp.T_SG_perc = T_SG_perc_curr
    new_hthp.calculate(calculate_thermo_only=True)

    return new_hthp


def opt_func(x, args=(1, 0.9)):

    print()
    print("{} - {}".format(x[0], args[1]))

    base_hthp = update_values(

        depth_in=depth,
        t_rock=t_rocks,
        t_amb=t_in_well,
        p_in=p_ref * x[0],
        T_SG_perc_curr=args[1]

    )

    min_value = 1 / (20 * base_hthp.m_dot_ratio) + args[0] * 3 / base_hthp.COP
    print("{} - {}   ->   {}".format(x[0], args[1], min_value))

    return abs(min_value)


# %%-------------------------------------   CHECK CONDITIONS                    -------------------------------------> #
x0 = [0.7, 0.9]
base_hthp = update_values(

    depth_in=depth,
    t_rock=t_rocks,
    t_amb=t_in_well,
    p_in=p_ref * x0[0],
    T_SG_perc_curr=x0[1]

)
print(base_hthp)
print(base_hthp.m_dot_ratio)
print(base_hthp.COP)


# %%-------------------------------------   RUN OPTIMIZATION                    -------------------------------------> #

x0 = [0.8]
x1s = np.linspace(0.1, 0.9, 5)
results = list()
opt_press = list()
for x1 in x1s:

    print(x1)
    res = minimize(opt_func, np.array(x0), args=[omega, x1], bounds=[(0.001, None)])
    results.append(res.fun)
    opt_press.append(res.x)

