# %%------------   IMPORT MODULES                     -----------------------------------------------------------> #
from main_code.simplified_well.simplified_well_subclasses import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot

# %% ------------ Number of Stages ----------------------------------------------------------->



#input condition
p_in = 12
t_in = 130

point_in=PlantThermoPoint(["CarbonDioxide"], [1])
point_in.set_variable("T",t_in)
point_in.set_variable("P",p_in)
h_in = point_in.get_variable("h")
s_in = point_in.get_variable("s")

#output condition
t_sat=28
point_sat=PlantThermoPoint(["CarbonDioxide"], [1])
point_sat.set_variable("T",t_sat)
point_sat.set_variable("Q",0)
p_out = point_sat.get_variable("p") #condenser pressure

point_out=PlantThermoPoint(["CarbonDioxide"], [1])
point_out.set_variable("P",p_out)
point_out.set_variable("s",s_in)
h_out = point_out.get_variable("h")

eta=0.8
h_out_r=h_in-(eta*(h_in-h_out))
point_out_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_out_r.set_variable("P",p_out)
point_out_r.set_variable("h",h_out_r)
s_out_r=point_out_r.get_variable("s")


delta_h_t=h_in-h_out
n=4 # initiation
delta_h_t_i=delta_h_t/n
P_t= 1000    # kW Fixing as desiered power output
m_dot= P_t/delta_h_t

point_tmp=PlantThermoPoint(["CarbonDioxide"], [1])
p_levels = list()

for i in range(n):

    point_tmp.set_variable("h", h_in-delta_h_t_i*(i+1))
    point_tmp.set_variable("s", s_in)
    p_levels.append(point_tmp.get_variable("P"))

def fi(rho,p,m=m_dot):
    Fi = m/((p*rho)**0.5)
    return Fi

def Y(Fi,B):
    Y= (1-(1/(B**2)))/Fi**2
    return Y

#1st stage
h_1_is=h_in-delta_h_t_i
point_1=PlantThermoPoint(["CarbonDioxide"], [1])
point_1.set_variable("h",h_1_is)
point_1.set_variable("s",s_in)
p_1 = point_1.get_variable("P")
rho_1 = point_1.get_variable("rho")

beta_1=p_in/p_1

h_1=h_in-(0.7*(h_in-h_1_is))

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
point_2.set_variable("s",s_in)
p_2 = point_2.get_variable("P")
rho_2 = point_2.get_variable("rho")
print("p_2",p_2)

h_2=h_1-(0.7*(h_1-h_2_is))
point_2_r=PlantThermoPoint(["CarbonDioxide"], [1])
point_2_r.set_variable("h",h_2)
point_2_r.set_variable("P",p_2)
s_2 = point_2_r.get_variable("s")
print("s_2",s_2)
beta_2=p_1/p_2
Fi_2= fi(rho_2,p_2,m_dot)
y_2=Y(Fi_2, beta_2)


#Third stage

h_3_is=h_2_is-delta_h_t_i
point_3=PlantThermoPoint(["CarbonDioxide"], [1])
point_3.set_variable("h",h_3_is)
point_3.set_variable("s",s_in)
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
point_4.set_variable("s",s_in)
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
print("y",y_1,y_2,y_3)
print("beta",beta_1,beta_2,beta_3)
print(p_in,"p_1=",p_1,"p_2",p_2,"p_3",p_3,"p_4",p_4,p_out)
print(p_levels)

# %% ------------ Design parameters ------------------------------------------------------------------->

Y_id=[0.4285451346283881, 0.34819140402298704, 0.2799601412216672, 0.22278180031350817]
beta_id=[1.1406090444111903, 1.1459273389990168, 1.1513816417219793, 1.1570025729163445]
Fi_id=[0.734751042618245, 0.8275795196845153, 0.9367591246345961, 1.06562533459034]




