# %%------------   IMPORT MODULES                     -----------------------------------------------------------> #
import math

import numpy as np

from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

# %% ------------ Design parameters ------------------------------------------------------------------->
Y_1_d=0.4285451346283881
Y_2_d=0.3481914040229870
Y_3_d=0.2799601412216672
Y_4_d=0.22278180031350817
Y=[0.4285451346283881,0.3481914040229870,0.2799601412216672,0.22278180031350817]
delta_h_1_d=11.917231959001676
delta_h_2_d=10.093400744332909
delta_h_3_d=8.47359080666627
delta_h_4_d=7.05088750743597
delta_h_i_d=[11.917231959001676,10.093400744332909,8.47359080666627,7.05088750743597]
Fi_1_D=0.7593907456076772
Fi_2_D=0.9014949852515453
Fi_3_D=1.0855694366899515
Fi_4_D=1.3259248989093555
Fi=[0.7593907456076772,0.9014949852515453,1.0855694366899515,1.3259248989093555]
T_amb=12
delta_T_appr=7
T_AC=T_amb+delta_T_appr
AC_point=PlantThermoPoint(["CarbonDioxide"], [1])
AC_point.set_variable("T",T_AC)
    if T_AC<30:
        AC_point.set_variable("Q", 0)
    else:
        AC_point.set_variable("rho", 700)
P_AC=AC_point.get_variable("P")

# %% ----------------stage 1------------------------------------------------------------------>
def Turbine_Expansion(m_r):
    T_in = 99.97
    P_in = 13.17

    Turbine_in_point = PlantThermoPoint(["CarbonDioxide"], [1])
    Turbine_in_point.set_variable("P", P_in)
    Turbine_in_point.set_variable("T", T_in)
    h_in = Turbine_in_point.get_variable("h")
    s_in = Turbine_in_point.get_variable("s")
    rho_in = Turbine_in_point.get_variable("rho")


    m_dot_D=50.12
    m_dot_OD=m_r*m_dot_D
    Fi_1 = m_dot_OD / np.sqrt(P_in * rho_in)
    beta_1 = 1 / (1 - (((Fi_1) ** 2) * Y_1_d)) ** 0.5

    P_1 = P_in / beta_1
    point_1 = PlantThermoPoint(["CarbonDioxide"], [1])
    point_1.set_variable("P", P_1)
    point_1.set_variable("s", s_in)
    T_1 = point_1.get_variable("T")
    h_1 = point_1.get_variable("h")
    rho_1 = point_1.get_variable("rho")

    Fi_2 = m_dot_OD / np.sqrt(P_1 * rho_1)
    beta_2 = 1 / np.sqrt(1 - (((Fi_2) ** 2) * Y_2_d))
    P_2 = P_1 / beta_2

    point_2 = PlantThermoPoint(["CarbonDioxide"], [1])
    point_2.set_variable("P", P_2)
    point_2.set_variable("s", s_in)
    h_2 = point_2.get_variable("h")
    T_2 = point_2.get_variable("T")
    rho_2 = point_2.get_variable("rho")

    Fi_3 = m_dot_OD / np.sqrt(P_2 * rho_2)
    beta_3 = 1 / (1 - (((Fi_3) ** 2) * Y_3_d)) ** 0.5
    P_3 = P_2 / beta_3
    point_3 = PlantThermoPoint(["CarbonDioxide"], [1])
    point_3.set_variable("P", P_3)
    point_3.set_variable("s", s_in)
    h_3 = point_3.get_variable("h")
    T_3 = point_3.get_variable("T")
    rho_3 = point_3.get_variable("rho")

    Fi_4 = m_dot_OD / np.sqrt(P_3 * rho_3)
    beta_4 = 1 / (1 - (((Fi_4) ** 2) * Y_4_d)) ** 0.5
    P_4 = P_3 / beta_4

    point_4 = PlantThermoPoint(["CarbonDioxide"], [1])
    point_4.set_variable("h", P_4)
    point_4.set_variable("s", s_in)
    h_4 = point_4.get_variable("h")
    T_4 = point_4.get_variable("T")
    rho_4 = point_4.get_variable("rho")

    beta_OD = beta_1 * beta_2 * beta_3 * beta_4
    #print(beta_OD, m_dot_OD)
    #print(P_AC)
    #print(P_4)
    return P_4
#print(Turbine_Expansion(0.9))
m_r=np.linspace(0.1,2.0,20)     # m vs P
for i in m_r:
    Turbine_Expansion(i)
    print('{}-{}-{}' .format(Turbine_Expansion(i), P_AC,i))
# %% ----------------iter------------------------------------------------------------------>


def Bisection(m_r_L,m_r_R):
    m_r_M = (m_r_R + m_r_L) / 2
    P_dif_M=Turbine_Expansion(m_r_M) - P_AC
    P_dif_L=Turbine_Expansion(m_r_L) - P_AC
    print('{}-{}'.format(P_dif_L,P_dif_M))

    iteration_counter=1

    if P_dif_M>0:
        if P_dif_L * P_dif_M > 0:
            m_r_L=m_r_M

        else:
            m_r_R=m_r_R+m_r_M

    else:
        if P_dif_L * P_dif_M > 0:   #both p_dif are negative, m_r_l should move towards left
            m_r_L=m_r_L-m_r_M

        else:
            m_r_R=m_r_M


    return P_dif_M, m_r_M

print(Bisection(0.8,0.9))





 # %% ----------------iter-
m_r=np.linspace(0.1,2.0,20)
print(len(m_r),m_r[0],m_r[19])

for a,b in m_r:
    a=float(m_r[0])
    b=float(m_r[19])
    Bisection(a,b)
    #print("{} - {} ".format(Bisection(i),Turbine_Expansion(i)))


#m_r=np.linspace(0.85,0.98,20)      prefered range




















