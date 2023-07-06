import math
from numpy import linspace
import numpy as np

from main_code.simplified_well.simplified_well_subclasses import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

# %% ------------ Design parameters ------------------------------------------------------------------->


delta_h_i_d=[11.917231959001676,10.093400744332909,8.47359080666627,7.05088750743597]


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
m_dot_D=50.12
def turbine_expansion(b_s):

    m_dot = b_s * m_dot_D



    T_in=100.5
    P_in=15.4

    tmp_point = PlantThermoPoint(["CarbonDioxide"], [1])
    tmp_point.set_variable("P",P_in)
    tmp_point.set_variable("T",T_in)
    h_in=tmp_point.get_variable("h")
    s_in=tmp_point.get_variable("s")
    rho_in=tmp_point.get_variable("rho")
    n_stages = 4
    points = list()             #turbine inlet / ist stage in
    points.append(tmp_point)
    #Y_ds = [Y_1_d, Y_2_d, Y_3_d, Y_4_d]
    Y = [0.4285451346283881, 0.3481914040229870, 0.2799601412216672, 0.22278180031350817]
    #Fi = [0.7593907456076772, 0.9014949852515453, 1.0855694366899515, 1.3259248989093555]

    for n in range(n_stages):
        P_st[n] = points[-1].get_variable("P")
        P_st.append(P_s t[n])

        rho_n = points[-1].get_variable("rho")
        rho_st.append(rho_n)

        s_n = points[-1].get_variable("s")
        s_st.append(s_n)

        Fi[n] = m_dot / np.sqrt(P_st[n] * rho_n)
        Fi.append(Fi[n])
        beta=list()
        beta_n = 1 / (1 - ((Fi_n ** 2) * Y[n])) ** 0.5
        beta.append(beta_n)
        #beta_T=(beta[0])*(beta[1])*(beta[2])*(beta[3])
        P_st[n] = float(P_st[n-1] )/ float(beta[n])

        tmp_point = tmp_point.duplicate()
        tmp_point.set_variable("P", P[4])
        tmp_point.set_variable("s", s_in)
        h_iso = tmp_point.get_variable("h")
        dh_iso = points[-1].get_variable("h") - h_iso
        print("{} - {}".format(P_out, beta_n))
        points.append(tmp_point)


    return points


# %%
turbine_expansion(1)


# %% ----------------iter------------------------------------------------------------------>


m_r_high = 2.5  #giving a range
m_r_low = 1
b_s = (m_r_high + m_r_low)/2
n_e=20
m_r=linspace(m_r_low,m_r_high,n_e)      #creating a range
print(m_r)


new_points = turbine_expansion(b_s)         #this returns a list of points
if math.isnan(points[-1].get_variable("P")):

    m_r_high = m_r

else:

    if AC_point.get_variable("P") - points[-1].get_variable("P") < 0:

        m_r_low = m_r

    else:

        m_r_high = m_r

print(AC_point.get_variable("P") - points[-1].get_variable("P"))




