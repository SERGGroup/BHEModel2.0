# %%------------   IMPORT MODULES                     -----------------------------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE
from tabulate import tabulate

# %%------------   BH input condition                   -----------------------------------------------------------> #

#let's assume a 0.7Mpa of pump pressure rise
delta_pump_p=0.7


#Pump input condition
T_amb = 15 # [°C]
dT_appr = 7  # [°C]
T_c_out=T_amb+dT_appr       #condenser Temperature


CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_c_out)


        if T_c_out < 30:

            CO2_input.set_variable("Q", 0.)
            p_c = CO2_input.get_variable("P")  # condenser pressure or turbine outlet pressure

        else:

            CO2_input.set_variable("rho", 700)
            p_c = CO2_input.get_variable("P")
        p_c = CO2_input.get_variable("P")
        print("condenser pressure",p_c)




CO2_input.set_variable("T", T_c_out)
p_in_BH=p_c+delta_pump_p
CO2_input.set_variable("P", p_in_BH)
h_in_BH = CO2_input.get_variable("h")

print("BH input pressure", p_in_BH)
print(CO2_input.get_unit("D"))





# %%------------     BH  Output condition               -----------------------------------------------------------> #

dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir

bhe_inside = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=dz_well, T_rocks=T_rock, use_rk=True

)


bhe_inside.update()
print(bhe_inside)

output_condition = bhe_inside.points[-1]
p_out_BH = output_condition.get_variable("P")
t_out_BH = output_condition.get_variable("T")
h_out_BH = output_condition.get_variable("h")
s_out_BH = output_condition.get_variable("s")
unit_p = output_condition.get_unit("P")
print(output_condition)

# %%------------ Turbine Expansion                      -----------------------------------------------------------> #

P_t_out=p_c

Turbine_input_condition=bhe_inside.points[-1]
p_t_in = Turbine_input_condition.get_variable("P")
t_t_in = Turbine_input_condition.get_variable("T")
h_t_in = Turbine_input_condition.get_variable("h")
s_t_in = Turbine_input_condition.get_variable("s")

Turbine_output_condition_iso = PlantThermoPoint(["CarbonDioxide"], [1])
s_t_out=Turbine_output_condition_iso.set_variable("s",s_t_in)
Turbine_output_condition_iso.set_variable("p", P_t_out)

h_out_iso = Turbine_output_condition_iso.get_variable("h")
t_t_out_iso= Turbine_output_condition_iso.get_variable("t")
print("turbine_out_t_iso:",t_t_out_iso)

eta=0.8

h_t_out=h_t_in-(eta*(h_t_in-h_out_iso))
print("h_t_out",h_t_out)
print("h_t_in",h_t_in)
delta_h_t=h_t_in-h_t_out
P_t= 1000                               # kW Fixing as desiered power output
m_dot= P_t/delta_h_t                    #kg/s
print("mass flow rate:Design",m_dot,"kg/s")

Turbine_output_condition_real = PlantThermoPoint(["CarbonDioxide"], [1])
Turbine_output_condition_real.set_variable("h",h_t_out)
Turbine_output_condition_real.set_variable("P", P_t_out)

s_t_out= Turbine_output_condition_real.get_variable("s")
t_t_out= Turbine_output_condition_real.get_variable("t")
print("turbine outlet temperature", t_t_out, "c")
beta=p_t_in/P_t_out
print("beta:",beta)

# %%------------ condenser                      -----------------------------------------------------------> #
T_c_in=t_t_out
condenser_output_condition= PlantThermoPoint(["CarbonDioxide"], [1])
condenser_output_condition.set_variable("Q",0)
condenser_output_condition.set_variable("P", p_c)

h_c_out= condenser_output_condition.get_variable("h")

print("h_c_out:",h_c_out)
Q_out=h_t_out-h_c_out
print("Q_out:",Q_out)

# %%------------   cycle efficiency neglecting pump wor ----------------------------------------------------------->
eta_cycle=100*(h_t_in-h_t_out)/(h_out_BH-h_in_BH)
print("eta_cycle:",eta_cycle)

# %%------------ Off design Analisys----------------------------------------------------------->
def massflow_coefficient(T,P,m=m_dot):
    Fi = m* (T**0.5)/P
    return Fi

def Y(P,Fi,B):
    Y= (P** 2 - B** 2) / ((P ** 2) * (Fi ** 2))
    return Y


Fi_d=massflow_coefficient(t_t_in,p_t_in,m_dot)
print("flow coefficient", Fi_d)

Y_id=Y(p_t_in,Fi_d,p_c)
print("Y_id",Y_id)

# %% ------------ Number of Stages ----------------------------------------------------------->
n=3 #just
delta_h_t_i=delta_h_t/n
P_3=p_t_in
h_3=h_t_in
s_3=s_t_in

h_2=h_3-delta_h_t_i
s_2
point_2=PlantThermoPoint(["CarbonDioxide"], [1])
point_2.set_variable("h",h_2)