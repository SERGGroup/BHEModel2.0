# %% INIT
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.off_design_model import TurbineOD, TurbineODPlotter
import matplotlib.pyplot as plt
import numpy as np

# %% DESIGN
from main_code.simplified_well.simplified_well_subclasses import SimplifiedBHE

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

# %% OFF-DESIGN
T_amb = 15 # [°C]
dT_appr = 7  # [°C]

beta_OD_list=list()
Fi_OD_list=list()
Power_OD_list=list()
m_dot_list=list()
T_a=list()
for T in range(0,24):

    T_amb=T
    T_a.append(T)
    turbine_in, turbine_out = evaluate_well()
    turbine.update_input_output_points(turbine_in, turbine_out)
    turbine.update_off_design_flow_rate()

    m_dot_od = turbine.input_point.m_dot
    m_dot_list.append(m_dot_od)

    h_in_od=turbine_in.get_variable("h")
    h_out_od=turbine_out.get_variable("h")
    rho_in=turbine_in.get_variable("rho")

    P_in=turbine_in.get_variable("P")
    P_out=turbine_out.get_variable("P")
    b=P_in/P_out
    Pow=(h_in_od-h_out_od)*m_dot_od
    Power_OD_list.append(Pow)
    beta_OD_list.append(b)
    fi=m_dot_od/(np.sqrt(P_in*rho_in))
    Fi_OD_list.append(fi)

# %% PLOT
plt.plot(T_a,Fi_OD_list)
plt.xlabel('T_amb/°c')
plt.ylabel('Mass Flow')
plt.show()
# %% PLOT
plt.plot(T_a,m_dot_list,'r')
plt.xlabel('T_amb/°c')
plt.ylabel('m_dot')
plt.show()
# %% PLOT
plt.plot(T_a,Power_OD_list,'g')
plt.xlabel('T_amb/°c')
plt.ylabel('Power/W')
plt.show()
# %% PLOT
plt.plot(T_a,beta_OD_list,'k')
plt.xlabel('T_amb/°c')
plt.ylabel('beta')
plt.show()