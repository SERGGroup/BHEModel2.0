# %%------------   IMPORT MODULES                     -----------------------------------------------------------> #
import math

import numpy as np

from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE

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
T_in=99.97
P_in=13.17

Turbine_in_point=PlantThermoPoint(["CarbonDioxide"], [1])
Turbine_in_point.set_variable("P",P_in)
Turbine_in_point.set_variable("T",T_in)
h_in=Turbine_in_point.get_variable("h")
s_in=Turbine_in_point.get_variable("s")
rho_in=Turbine_in_point.get_variable("rho")


m_dot_D=50.12
m_r=0.9
m_dot_OD=m_r*m_dot_D
Fi_1=m_dot_OD/np.sqrt(P_in*rho_in)
beta_1=1/(1-(((Fi_1)**2)*Y_1_d))**0.5

P_1=P_in/beta_1
point_1=PlantThermoPoint(["CarbonDioxide"], [1])
point_1.set_variable("P",P_1)
point_1.set_variable("s",s_in)
T_1=point_1.get_variable("T")
h_1=point_1.get_variable("h")
rho_1=point_1.get_variable("rho")



# %% ----------------stage 2------------------------------------------------------------------>
Fi_2=m_dot_OD/np.sqrt(P_1*rho_1)
beta_2=1/np.sqrt(1-(((Fi_2)**2)*Y_2_d))
P_2=P_1/beta_2

point_2=PlantThermoPoint(["CarbonDioxide"], [1])
point_2.set_variable("P",P_2)
point_2.set_variable("s",s_in)
h_2=point_2.get_variable("h")
T_2=point_2.get_variable("T")
rho_2=point_2.get_variable("rho")


# %% ----------------stage 3------------------------------------------------------------------>
Fi_3=m_dot_OD/np.sqrt(P_2*rho_2)
beta_3=1/(1-(((Fi_3)**2)*Y_3_d))**0.5
P_3=P_2/beta_3
point_3=PlantThermoPoint(["CarbonDioxide"], [1])
point_3.set_variable("P",P_3)
point_3.set_variable("s",s_in)
h_3=point_3.get_variable("h")
T_3=point_3.get_variable("T")
rho_3=point_3.get_variable("rho")


# %% ----------------stage 4------------------------------------------------------------------>
Fi_4=m_dot_OD/np.sqrt(P_3*rho_3)
beta_4=1/(1-(((Fi_4)**2)*Y_4_d))**0.5
P_4=P_3/beta_4

point_4=PlantThermoPoint(["CarbonDioxide"], [1])
point_4.set_variable("h",P_4)
point_4.set_variable("s",s_in)
h_4=point_4.get_variable("h")
T_4=point_4.get_variable("T")
rho_4=point_4.get_variable("rho")


beta_OD=beta_1*beta_2*beta_3*beta_4
print(beta_OD,m_dot_OD)
print(P_AC)
print(P_4)

# %% ----------------iter------------------------------------------------------------------>

def turbine_expansion(m_r):

    m_dot = m_r * m_dot_D

    n_stages = 4

    T_in=99.97
    P_in=13.17

    tmp_point = PlantThermoPoint(["CarbonDioxide"], [1])
    tmp_point.set_variable("P",P_in)
    tmp_point.set_variable("T",T_in)

    points = list()
    points.append(tmp_point)
    Y_ds = [Y_1_d, Y_2_d, Y_3_d, Y_4_d]
    h_d=[delta_h_1_d,delta_h_2_d,delta_h_3_d,delta_h_4_d]
    for n in range(n_stages):

        P_in = points[-1].get_variable("P")
        rho_in = points[-1].get_variable("rho")
        h_in = points[-1].get_variable("h")
        s_in = points[-1].get_variable("s")

        Fi_n = m_dot / np.sqrt(P_in * rho_in)

        beta_n = 1 / (1 - ((Fi_n ** 2) * Y[n])) ** 0.5
        #a_n = np.log(delta_h_n / h_[n])
        P_out = P_in / beta_n

        tmp_point = tmp_point.duplicate()
        tmp_point.set_variable("P", P_out)
        tmp_point.set_variable("s", s_in)
        h_iso = tmp_point.get_variable("h")
        dh_iso = points[-1].get_variable("h") - h_iso
        print(P_out)

        print("{} - {}".format(P_out, beta_n))

        points.append(tmp_point)

    return points

# %% ----------------iter------------------------------------------------------------------>
import math

m_dot_D=50.12
m_r_high = 5
m_r_low = 0


m_r =np.linspace(m_r_low,m_r_high,20)
for i in m_r:
    print(turbine_expansion(i))
#points = turbine_expansion(m_r)

#if math.isnan(points[-1].get_variable("P")):

    #m_r_high = m_r

#else:

    #if AC_point.get_variable("P") - points[-1].get_variable("P") < 0:

        #m_r_low = m_r

    #else:

        #m_r_high = m_r

#print(AC_point.get_variable("P") - points[-1].get_variable("P"))




# %% ----------------efficincy calculation------------------------------------------------------------------>
delta_h_1=h_in-h_1
delta_h_2=h_1-h_2
delta_h_3=h_2-h_3
delta_h_4=h_3-h_4

a_1=np.log(delta_h_1/delta_h_1_d)
a_2=np.log(delta_h_2/delta_h_2_d)
a_3=np.log(delta_h_3/delta_h_3_d)
a_4=np.log(delta_h_4/delta_h_4_d)

eta_1=0.7*(10**((-0.00817*(a_1**3))-(0.03181*(a_1**2))+0.0019*a_1))
eta_2=0.7*(10**((-0.00817*(a_2**3))-(0.03181*(a_2**2))+0.0019*a_2))
eta_1=0.7*(10**((-0.00817*(a_3**3))-(0.03181*(a_3**2))+0.0019*a_3))
eta_1=0.7*(10**((-0.00817*(a_4**3))-(0.03181*(a_4**2))+0.0019*a_4))
print(eta_1)
 