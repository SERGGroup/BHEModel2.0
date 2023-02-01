# %% ------------
import math
from numpy import linspace
import numpy as np

from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE
# %% ------------Design parameters ------------------------------------------------------------------->
m_dot_D=43.85583
delta_h_i_d=[4.987933577345359, 6.484313650548927, 6.933227672510043 ,7.067901879098372]
Y_ds=[0.4525192696952781, 0.3582619774594879, 0.27740098563753757, 0.21005967202839923]
#OD condition
T_amb=29
delta_T_appr=7
T_AC=T_amb+delta_T_appr
AC_point=PlantThermoPoint(["CarbonDioxide"], [1])
AC_point.set_variable("T",T_AC)
    if T_AC<30:
        AC_point.set_variable("Q", 0)
    else:
        AC_point.set_variable("rho", 700)
P_AC=AC_point.get_variable("P")
print("P_AC",P_AC)

# %% ----------------Turbine expansion------------------------------------------------------------------>

#n_e=20
#m_r=linspace(m_r_low,m_r_high,n_e)      #slicing method

def turbine_expansion(m_r):

    m_r_high = 2.5  # giving a range
    m_r_low = 1
    m_r = (m_r_high + m_r_low) / 2
    m_dot = m_r * m_dot_D

    P = [13.17]
    T = [99.97]
    eta_des = 0.7

    tmp_point = PlantThermoPoint(["CarbonDioxide"], [1])
    tmp_point.set_variable("P", P[0])
    tmp_point.set_variable("T", T[0])

    rho = [tmp_point.get_variable("rho")]
    h = [tmp_point.get_variable("h")]
    s = [tmp_point.get_variable("s")]
    points = [tmp_point]
    h_iso = list()

    Fi = list()
    beta = list()

    for n in range(len(Y_ds)):  # 0 1 2 3

        P_in = P[-1]
        rho_in = rho[-1]
        s_in = s[-1]
        h_in = h[-1]

        Fi_curr = m_dot / (np.sqrt(P_in * rho_in))
        beta_curr = 1 / (1 - ((Fi_curr ** 2) * Y_ds[n])) ** 0.5
        P_out = P_in/beta_curr

        tmp_point = tmp_point.duplicate()
        tmp_point.set_variable("P", P_out)
        tmp_point.set_variable("s", s_in)

        h_iso_curr = tmp_point.get_variable("h")
        delta_h_iso = h_in - h_iso_curr
        a = np.log(delta_h_iso / delta_h_i_d[n])
        eta_curr = eta_des * (10 ** ((-0.00817 * (a ** 3)) - (0.03181 * (a ** 2)) + 0.0019 * a))

        h_out = h_in - eta_curr * delta_h_iso

        tmp_point.set_variable("P", P_out)
        tmp_point.set_variable("h", h_out)

        Fi.append(Fi_curr)
        beta.append(beta_curr)
        P.append(P_out)
        rho.append(tmp_point.get_variable("rho"))
        h_iso.append(h_iso_curr)
        points.append(tmp_point)

    return points

print(turbine_expansion((0.5)))

 # %% ----------------efficincy calculation------------------------------------------------------------------>
for i in m_r:
    turbine_expansion(i)
    print('{}-{}-{}' .format(turbine_expansion(i), P_AC,i))

def Bisection(m_r_L,m_r_R):
    m_r_M = (m_r_R + m_r_L) / 2
    P_dif_M=turbine_expansion(m_r_M) - P_AC
    P_dif_L=turbine_expansion(m_r_L) - P_AC
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


print(P_dif_M, m_r_M)

print(Bisection(0.6,0.7))


