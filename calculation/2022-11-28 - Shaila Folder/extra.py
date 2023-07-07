#%%------------   IMPORT MODULES                     -----------------------------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot

# %%------------   BH input condition                   -----------------------------------------------------------> #

# let's assume a 0.7Mpa of pump pressure rise
delta_pump_p = 0.7


# Pump input condition
T_amb = 15 # [°C]
dT_appr = 7  # [°C]
T_c_out=T_amb+dT_appr       # condenser Temperature


CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_c_out)


if T_c_out < 30:        # Transcritical

    CO2_input.set_variable("Q", 0.)
    p_c = CO2_input.get_variable("P")  # condenser pressure or turbine outlet pressure
else:                   # Supercritical

    CO2_input.set_variable("rho", 700)
    p_c = CO2_input.get_variable("P")
    ip_c = CO2_input.get_variable("P")


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
    dz_well=dz_well, t_rocks=T_rock, use_rk=True

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

eta=0.7

h_t_out=h_t_in-(eta*(h_t_in-h_out_iso))
print("h_t_out",h_t_out)
print("h_t_in",h_t_in)
delta_h_t=h_t_in-h_t_out
delta_h_t_iso=h_t_in-h_out_iso
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

# %%------------ condenser   / Air cooler                   -----------------------------------------------------------> #
T_c_in=t_t_out

condenser_output_condition= PlantThermoPoint(["CarbonDioxide"], [1])
condenser_output_condition.set_variable("Q",0)
condenser_output_condition.set_variable("P", p_c)
h_c_out = condenser_output_condition.get_variable("h")



print("h_c_out:",h_c_out)
Q_out=h_t_out-h_c_out
print("Q_out:",Q_out)

# %%------------   cycle efficiency neglecting pump wor ----------------------------------------------------------->
eta_cycle=100*(h_t_in-h_t_out)/(h_out_BH-h_in_BH)
print("eta_cycle:",eta_cycle)

# %%------------ Off design Analisys----------------------------------------------------------->
def fi(rho,p,m=m_dot):
    Fi = m/((p*rho)**0.5)
    return Fi

def Y(Fi,B):
    Y= (1-(1/(B**2)))/Fi**2
    return Y





# %% ------------ Number of Stages ----------------------------------------------------------->
#Assuming constant enthalpy drop accros each stage


n=4 #just
delta_h_t_i=delta_h_t_iso/n
point_tmp=PlantThermoPoint(["CarbonDioxide"], [1])
p_levels = list()

for i in range(n):

    point_tmp.set_variable("h", h_t_in-delta_h_t_i*(i+1))
    point_tmp.set_variable("s", s_t_in)
    p_levels.append(point_tmp.get_variable("P"))

#1st stage
h_1_is=h_t_in-delta_h_t_i

point_1=PlantThermoPoint(["CarbonDioxide"], [1])
point_1.set_variable("h",h_1_is)
point_1.set_variable("s",s_t_in)
p_1 = point_1.get_variable("P")
rho_1 = point_1.get_variable("rho")

beta_1=p_t_in/p_1

h_1=h_t_in-(0.7*(h_t_in-h_1_is))

point_1_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_1_r.set_variable("h",h_1)
point_1_r.set_variable("P",p_1)
s_1 = point_1_r.get_variable("s")
Fi_1= fi(rho_1,p_1,m_dot)
y_1=Y(Fi_1, beta_1)

#second stage

h_2_is=h_1_is-delta_h_t_i
point_2=PlantThermoPoint(["CarbonDioxide"], [1])
point_2.set_variable("h",h_2_is)
point_2.set_variable("s",s_t_in)
p_2 = point_2.get_variable("P")
rho_2 = point_2.get_variable("rho")

h_2=h_1-(0.7*(h_1-h_2_is))
point_2_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_2_r.set_variable("h",h_2)
point_2_r.set_variable("P",p_2)
s_2 = point_2_r.get_variable("s")
beta_2=p_1/p_2
Fi_2= fi(rho_2,p_2,m_dot)
y_2=Y(Fi_2, beta_2)


#Third stage

h_3_is=h_2_is-delta_h_t_i
point_3=PlantThermoPoint(["CarbonDioxide"], [1])
point_3.set_variable("h",h_3_is)
point_3.set_variable("s",s_t_in)
p_3 = point_3.get_variable("P")
rho_3 = point_3.get_variable("rho")

h_3=h_2-(0.7*(h_2-h_3_is))
point_3_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_3_r.set_variable("h",h_3)
point_3_r.set_variable("P",p_3)
s_3 = point_3_r.get_variable("s")
beta_3=p_2/p_3
Fi_3= fi(rho_3,p_3,m_dot)
y_3=Y(Fi_3, beta_3)

#4th stage
h_4_is=h_3_is-delta_h_t_i
point_4=PlantThermoPoint(["CarbonDioxide"], [1])
point_4.set_variable("h",h_4_is)
point_4.set_variable("s",s_t_in)
p_4 = point_4.get_variable("P")
rho_3 = point_3.get_variable("rho")


h_4=h_3-(0.7*(h_3-h_4_is))
point_4_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_4_r.set_variable("h",h_4)
point_4_r.set_variable("P",p_4)
s_4 = point_4_r.get_variable("s")
rho_4 = point_4.get_variable("rho")
beta_4=p_3/p_4
Fi_4= fi(rho_4,p_4,m_dot)
y_4=Y(Fi_4, beta_4)
print("y",y_1,y_2,y_3,y_4)
print("beta",beta_1,beta_2,beta_3,beta_4)
print("p_t_in=",p_t_in,"p_1=",p_1,"p_2",p_2,"p_3",p_3,"p_4", p_4,"p_c",p_c)
print(p_levels)
print("Fi", Fi_1,Fi_2,Fi_3,Fi_4)


