# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint


# %%------------   INIT INPUT CONDITION                   -----------------------------------------------------------> #
T_amb = 22 # [°C]
dT_appr = 7  # [°C]
Y_id = 5213221

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_amb + dT_appr)

if T_amb + dT_appr < 31.1:

    CO2_input.set_variable("Q", 0.)
    p_in = CO2_input.get_variable("P")

else:

    CO2_input.set_variable("rho", 700)
    p_in = CO2_input.get_variable("P")

CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("P", p_in*1.0001)

print("BH input pressure", p_in)
print(CO2_input.get_unit("D"))

 # %%------------   INIT BHE WELL                          -----------------------------------------------------------> #
dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir

bhe_in = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=dz_well, t_rocks=T_rock, use_rk=True

)


bhe_in.update()
print(bhe_in)

output_condition = bhe_in.points[-1]
p_out = output_condition.get_variable("P")
t_out = output_condition.get_variable("T")
h_out = output_condition.get_variable("h")
s_out = output_condition.get_variable("s")
unit_p = output_condition.get_unit("P")
print(output_condition)

# %%------------   TURBINE EXPANSION                      -----------------------------------------------------------> #
delta_P_pump=0.7  #pump pressure rise
P_t_out=p_in-delta_P_pump
print(P_t_out)

Turbine_input_condition=bhe_in.points[-1]
p_t_in = Turbine_input_condition.get_variable("P")
t_t_in = Turbine_input_condition.get_variable("T")
h_t_in = Turbine_input_condition.get_variable("h")
s_t_in = Turbine_input_condition.get_variable("s")

Turbine_output_condition_iso = PlantThermoPoint(["CarbonDioxide"], [1])
s_t_out=Turbine_output_condition_iso.set_variable("s",s_t_in)
Turbine_output_condition_iso.set_variable("p", P_t_out)

h_out_iso = Turbine_output_condition_iso.get_variable("h")
t_t_out_iso= Turbine_output_condition_iso.get_variable("t")
print("t_out_t_iso:",t_t_out_iso)
eta=0.8
h_t_out=h_t_in-(eta*(h_t_in-h_out_iso))
print("h_t_out",h_t_out)
print("h_t_in",h_t_in)
delta_h_t=h_t_in-h_t_out
P_t= 1000  # kW
m_dot= P_t/delta_h_t   #kg/s
print("mass flow rate:Design",m_dot,"kg/s")

Turbine_output_condition_real = PlantThermoPoint(["CarbonDioxide"], [1])
Turbine_output_condition_real.set_variable("h",h_t_in)
Turbine_output_condition_real.set_variable("P", P_t_out)

s_t_out= Turbine_output_condition_real.get_variable("s")
t_t_out= Turbine_output_condition_real.get_variable("t")
print("turbine outlet temperature", t_t_out, "c")
beta=p_t_in/P_t_out
print("beta:",beta)

# %%------------   Condensser                      ----------------------------------------------------------->

p_c= p_in-delta_P_pump  #condenser pressure
delta_t_c=t_t_out-(T_amb+dT_appr)
condenser_output_condition= PlantThermoPoint(["CarbonDioxide"], [1])
condenser_output_condition.set_variable("Q",0)
condenser_output_condition.set_variable("T", T_amb + dT_appr)


h_c_out= condenser_output_condition.get_variable("h")
t_c_out= condenser_output_condition.get_variable("t")
print(t_c_out)
print(h_c_out)
Q_out=h_t_out-h_c_out
print("Q_out:",Q_out)
# %%------------   Pump work                      ----------------------------------------------------------->
Pump_p_in=p_c
Pump_W=h_out-h_c_out
print("Pump_W:",Pump_W)

# %%------------   Cycle efficiency                      ----------------------------------------------------------->
eta_des=0.7

# %%------------   evaluate off design                     ----------------------------------------------------------->

def calc_function(t_amb, opt_parm):

    # Set the input condition for the BHE well

    # evaluate the turbine off_desing


    # evauate and return the power