# %% INIT
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.off_design_model import TurbineOD, TurbineODPlotter
import matplotlib.pyplot as plt

# %% DESIGN
from main_code.simplified_well.simplified_well import SimplifiedBHE

T_amb = 15 # [°C]
dT_appr = 7  # [°C]
delta_pump_p = 0.1

dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir


def evaluate_well():

    T_c_out=T_amb+dT_appr       # condenser Temperature
    CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
    CO2_input.set_variable("T", T_c_out)
    CO2_input.set_variable("Q", 0.)

    p_c = CO2_input.get_variable("P")
    CO2_input.set_variable("T", T_c_out)
    p_in_BH = p_c + delta_pump_p
    CO2_input.set_variable("P", p_in_BH)

    bhe_inside = SimplifiedBHE(

        input_thermo_point=CO2_input,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

    )

    bhe_inside.update()

    return bhe_inside.points[-1].duplicate(), CO2_input.duplicate()

turbine_in, turbine_out = evaluate_well()
turbine = TurbineOD(turbine_in, turbine_out, n_stages=3)

# %% INPUT POINTS

conditions = list()
for T in range(-5, 20):

    T_amb = T
    turbine_in, turbine_out = evaluate_well()
    conditions.append([turbine_in, turbine_out])


# %% PLOT

fig = plt.figure(dpi=150)
fig.set_size_inches(20, 8)
turbine_plotter = TurbineODPlotter(turbine)
turbine_plotter.plot_characteristic_curves(fig, n_points=40, off_design_conditions=conditions)
plt.tight_layout(pad=5)
plt.show()
