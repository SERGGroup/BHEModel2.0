# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.power_plants.base_surface_power_plant import BaseCO2HTHP
from scipy.optimize import minimize
import numpy as np
import gc


# %%-------------------------------------   INITIALIZE CALCULATIONS             -------------------------------------> #
t_in_well = 10
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


def update_values(depth_in, t_amb, p_in):

    t_rocks = t_amb + grad * depth / 1000

    bhe_in.set_variable("T", t_amb)
    bhe_in.set_variable("P", p_in)

    bhe_well = SimplifiedBHE(

        bhe_in, depth_in,
        t_rocks=t_rocks, t_surf=t_amb

    )
    new_hthp = BaseCO2HTHP(

        bhe_well=bhe_well,
        p_steam=p_steam,
        p_max=50

    )
    new_hthp.prepare_calculation()
    new_hthp.calculate()

    return new_hthp


def opt_func(x, args=(1, 20)):

    base_hthp = update_values(

        depth_in=depth,
        t_amb=args[1],
        p_in=p_ref * x[0]

    )

    if base_hthp.is_physical:

        min_value = 1 / (20 * base_hthp.m_dot_ratio) + args[0] * 3 / base_hthp.COP
        print("{} - {} -> {}".format(x[0], args[1], min_value))

        return abs(min_value)

    return np.inf


# %%-------------------------------------   CHECK CONDITIONS                    -------------------------------------> #
x0 = [0.1]
base_hthp = update_values(

    depth_in=depth,
    t_amb=10,
    p_in=p_ref * x0[0],

)
print(base_hthp)
print(base_hthp.m_dot_ratio)
print(base_hthp.COP)
p_low = base_hthp.points[3].get_variable("P")
tmp_point = base_hthp.points[0].duplicate()
tmp_point.set_to_expansion_result(p_low, 0.75, base_hthp.points[2])
print(tmp_point.get_variable("T"))
print(base_hthp.is_p_max_limited)


# %%-------------------------------------   RUN OPTIMIZATION                    -------------------------------------> #
t_in_wells = np.linspace(-10, 30, 5)
results = list()
opt_press = list()
for t_in_well in t_in_wells:

    print(t_in_well)
    res_vap = minimize(opt_func, np.array([0.01]), args=[omega, t_in_well], bounds=[(0.001, 0.55)])
    results.append(res_vap.fun)
    opt_press.append(res_vap.x)

