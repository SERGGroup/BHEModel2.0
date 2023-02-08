# %%------------   IMPORT MODULES                     -----------------------------------------------------------> #
from main_code.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

# %% ------------ Design parameters ------------------------------------------------------------------->
Y_1_d=0.4285451346283881
Y_2_d=0.3481914040229870
Y_3_d=0.2799601412216672
Y_4_d=0.22278180031350817
delta_h_1_d=11.917231959001676
delta_h_2_d=10.093400744332909
delta_h_3_d=8.47359080666627
delta_h_4_d=7.05088750743597
m_dot=36.38  # required mass flow to produce 1MW at 22 degrees
#when ambient temperature changes the mass flowrate remain constant/ turbine inlet and outlet condition changes


T_amb=[0	,1	,2	,5	,8	,10	,12	,14,17	,20,22	,23	,24	,25	,26	,27,28,	29,	30]
delta_T_appr=7
T_c=list()
for element in T_amb:
    value = int(element)+delta_T_appr
    T_c.append(value)
print(T_c)
#Turbine inlet conditions/ well output condition
T_1=[99.98,99.98,99.98,99.98,99.97,99.97,99.97,99.97,99.97,99.98,99.99,100.01,100.1,100.16,100.23,100.31,100.4,100.5,100.61]
P_1=[13.08,13.09,13.1,13.13,13.15,13.17,13.17,13.17,13.14,13.04,12.87,12.69,14.26,14.48,14.71,14.93,15.17,15.4,15.63]



AC_pressure=list()
AC_point=PlantThermoPoint(["CarbonDioxide"], [1])
for element in T_c:
    a=int(element)
    AC_point.set_variable("T",a)
    if a<30:
        AC_point.set_variable("Q", 0)
    else:
        AC_point.set_variable("rho", 700)
    AC_pressure.append(AC_point.get_variable("P"))


# %% ----------------stage 1------------------------------------------------------------------>
Turbine_in_point=PlantThermoPoint(["CarbonDioxide"], [1])
Turbine_in_h=list()
Turbine_in_s=list()
Turbine_in_rho=list()
for element in AC_pressure:
    Turbine_in_point.set_variable("P",float(element))
    for element in T_c:
        Turbine_in_point.set_variable("T",int(element))
        Turbine_in_h.append(Turbine_in_point.get_variable("h"))
        Turbine_in_s.append(Turbine_in_point.get_variable("s"))
        Turbine_in_rho.append(Turbine_in_point.get_variable("rho"))

Fi_Turbine_in=list()

for element in Turbine_in_rho:
    x=float(element)
    for element in P_1:
        y=float(element)
        z=m_dot/((x*y)**0.5)
        Fi_Turbine_in.append(z)

beta_1=list()
for element in Fi_Turbine_in:
    b=element
    d=1/(1-(((b)**2)*Y_1_d))**0.5
    beta_1.append(d)


# %% ----------------stage 2------------------------------------------------------------------>
P_2=list()
for element in beta_1:
    e=float(element)
    for element in P_1:
        f=float(element)
        g=f/e
P_2.append(g)

h_2=list()      #isentropic
rho_2=list()
T_2=list()
point_1=PlantThermoPoint(["CarbonDioxide"], [1])
for element in P_2:
    j=float(element)
    point_1.set_variable("P",j )
    for element in Turbine_in_s:
        k=float(element)
        point_1.set_variable("s", k)
        l=point_1.get_variable("T")
        m=point_1.get_variable("h")
        n=point_1.get_variable("rho")
        h_2.append(m)
        T_2.append(l)
        rho_2.append(n)

delta_h_1=list()
for element in Turbine_in_h:
    x=element
    for elementin in h_2:
        y=element
        z=x-y
        delta_h_1.append(z)


# %% ----------------efficiency curve------------------------------------------------------------------>



