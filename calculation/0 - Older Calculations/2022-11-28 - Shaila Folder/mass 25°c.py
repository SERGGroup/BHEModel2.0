from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.power_plants.off_design_model.turbine_off_design import TurbineOD
from pathlib import Path
import pandas as pd
from pandas import *
import numpy as np
from numpy import linspace
from pandas import Series,DataFrame
import matplotlib.pyplot as plt
from tqdm import tqdm
import calendar
from datetime import datetime, date
from tabulate import tabulate
import pandas as pd

# %%
city=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Munich ERA.csv'
#city=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Florence ERA.csv'
#city=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Cheyenne ERA.csv'
#city=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Riyadh ERA.csv'
meteo=pd.read_csv(city,usecols=['time', 'temperature','relativehumidity','pressure'])
average_T=''
average_H=''
average_P=''
meteo['average_T']=average_T
meteo['average_H']=average_H
meteo['average_P']=average_P
meteo['time'] = pd.to_datetime(meteo['time'])
Monthly_df=meteo[meteo["time"].between
        ('2018-01-01T00:00', '2018-12-31T00:00')]
T_av=Monthly_df.groupby([Monthly_df['time'].dt.month,
                                    Monthly_df['time'].dt.hour]).temperature.mean()
RH_av=Monthly_df.groupby([Monthly_df['time'].dt.month,
                                    Monthly_df['time'].dt.hour]).relativehumidity.mean()
P_av=Monthly_df.groupby([Monthly_df['time'].dt.month,
                                    Monthly_df['time'].dt.hour]).pressure.mean()
print(P_av)
# %%


#powr = tqdm(desc="Calculation", total=len(T_av))
annual_T_list=list(T_av)
mon = linspace(0, 11, 12)
AC_pressure = list()
T_Cond_annual = list()
h_c_out = list()
BH_in_point = list()
BH_out_points = list()
BH_in_point_OD = list()
BH_out_point_OD = list()
h_in = list()  # turbine in
s_in = list()
p_in=list()
rho_in=list()
heat_sink = list()
heat_source = list()
cycle_eta = list()
delta_T_appr = 10
Turbine_out_r_points = list()   #turbine out
eta = 0.7
T_out = list
h_out_iso = list()
h_out_r = list()
work = list()
Turbine_Power = list()
pressure_ratio=list()
fi_in=list()
c_in=list()
Ma_in=list()
vol_in=list()
point = PlantThermoPoint(["Carbon Dioxide"], [1])
for i in annual_T_list:
    T_mean = float(i) + delta_T_appr
    T_Cond_annual.append(T_mean)


for i in T_Cond_annual:
    tmp_point = point.duplicate()
    tmp_point.set_variable("T", i)
    if i < 30:
        tmp_point.set_variable("Q", 0)
    else:
        tmp_point.set_variable("rho", 700)
    AC_pressure.append(tmp_point.get_variable("P"))
    BH_in_point.append(tmp_point)


# %%

BH_points = PlantThermoPoint(["Carbon Dioxide"], [1])

dz_well = 1500  # [m] Depth of the reservoir                    #BH output T and P
T_rock = 125  # [°C] Temperature of the reservoir
for element in BH_in_point:
    bhe_in = SimplifiedBHE(

        input_thermo_point=element,
        dz_well=dz_well, t_rocks=T_rock, use_rk=True

    )

    bhe_in.update()
    output_condition = bhe_in.points[-1]
    h_in.append(output_condition.get_variable("h"))
    s_in.append(output_condition.get_variable("s"))
    p_in.append(output_condition.get_variable("p"))
    rho_in.append(output_condition.get_variable("rho"))
    c_in.append(output_condition.get_variable("c"))

    bhe_in.points[-1].copy_state_to(BH_points)
    BH_out_points.append(BH_points)

TO_point = PlantThermoPoint(["Carbon Dioxide"], [1])
for i in range(len(AC_pressure)):
    tmp_point = TO_point.duplicate()
    tmp_point.set_variable("s", s_in[i])
    tmp_point.set_variable("P", AC_pressure[i])
    h_out_iso.append(tmp_point.get_variable("h"))


# %%
ac_point = PlantThermoPoint(["Carbon Dioxide"], [1])
for y in T_Cond_annual:
    tmp_point = ac_point.duplicate()
    tmp_point.set_variable("T", y)
    if y < 30:
        tmp_point.set_variable("Q", 0)
    else:
        tmp_point.set_variable("rho", 700)

    h_c_out.append(tmp_point.get_variable("h"))




# %% design
turb_des_out = PlantThermoPoint(["Carbon Dioxide"], [1])

t_des= 25
turb_des_out.set_variable("T", t_des)

if t_des < 30:
    turb_des_out.set_variable("Q", 0)

else:
    turb_des_out.set_variable("rho", 700)

P_cond_d=turb_des_out.get_variable("P")



# %%
turb_des_in = PlantThermoPoint(["Carbon Dioxide"], [1])
dz_well = 1500
T_rock = 125

bhe_in = SimplifiedBHE(

    input_thermo_point=turb_des_out,
    dz_well=dz_well, t_rocks=T_rock, use_rk=True)

bhe_in.update()
bhe_in.points[-1].copy_state_to(turb_des_in)
turb_des_in.m_dot = 76.70
for i in range(len(AC_pressure)):
    pressure_ratio.append(p_in[i]/AC_pressure[i])

# %%

Turbine_Power_OD = list()
Turbine_m_dot_OD = list()
Turbine_eta_OD = list()

turbine_d=list()
turbine = TurbineOD(input_point=turb_des_in, output_point=turb_des_out,n_stages=3)
turbine.evaluate_design()



dh_iso_1=list()
fi=list()
for j in range(len(BH_in_point)):
    k = BH_in_point[j]
    l = BH_out_points[j]
    turbine.update_input_output_points(l, k)
    turbine.update_off_design_flow_rate()
    Turbine_Power_OD.append(turbine.power)
    Turbine_m_dot_OD.append(turbine.input_point.m_dot)
    Turbine_eta_OD.append(turbine.eta_iso)

for j in range(len(Turbine_eta_OD)):       #this is wrong
    e = h_in[j]
    f = h_out_iso[j]
    x = e - (Turbine_eta_OD[j] * (e - f))
    h_out_r.append(x)
    work.append(e-x)

# %%
for q in range(len(h_out_r)):
    r = h_out_r[q]
    l = h_c_out[q]
    o = h_in[q]
    n = r - l
    u = o - l
    heat_sink.append(h_out_r[q]-h_c_out[q])
    heat_source.append(h_in[q]-h_c_out[q])

# %%
for w in range(len(heat_source)):
    x = heat_source[w]
    y = work[w]
    eta = (y / x) * 100
    cycle_eta.append(eta)

for i in range(len(h_in)):
    fi_in.append(Turbine_m_dot_OD[i]/((rho_in[i]*p_in[i])**0.5))




# %%
heat_sink_tot=list()
for i in range(len(heat_sink)):
    heat_sink_tot.append(Turbine_m_dot_OD[i]*heat_sink[i])

# %%----- mass flowrate vs t_amb


fig_8,ax3 = plt.subplots()
ax3.plot(annual_T_list,Turbine_m_dot_OD,color= 'r')
ax3.plot(annual_T_list,[76.70]*288,label="m_des")
plt.ylabel("Mass flow rate/ kg/s")
plt.xlabel("Ambient temperature")
plt.title("Mass flow rate variation with Ambient Temperature (T_des=25°c)")
#plt.title("Mass flow rate variation with Ambient Temperature,Florence(T_des=1°c)")
#plt.title("Mass flow rate variation with Ambient Temperature,Cheyenne(T_des=1°c)")
#plt.title("Mass flow rate variation with Ambient Temperature,Riyadh(T_des=1°c)")
plt.ylim(10,140)
plt.legend()
plt.grid()
plt.show()

# %%----- Turb inlet/outlet enthalpy--------------
a=list(h_in)
b=list(h_out_r)
fig_2,ax=plt.subplots()
ax.plot(annual_T_list,a,label="in")
ax.plot(annual_T_list,b,label="out")
ax.set_ylabel("Enthalpy (kJ/kg)")
ax.set_xlabel("Ambient Temperature (°C)")
ax.set_ylim(420,500)
ax.set_title("Turbine inlet/outlet enthalpy variation with T_amb(T_des=25°C)" )
ax.legend()
plt.show()
# %%----- BH inlet/outlet pressure, pressure ratio
fig_3,ax1=plt.subplots()

ax1.plot(annual_T_list,AC_pressure,label="P_out")
ax1.plot(annual_T_list,p_in,label="P_in")
ax1.set_xlabel("Ambient Temperature/ °c")
ax1.set_ylabel("Pressure/ MPa")
ax1.set_ylim(0,18)
ax1.legend(loc='lower left')
ax2=ax1.twinx()
ax2.plot(annual_T_list,pressure_ratio,label="Pressure ratio",color="g")
ax2.set_ylim(0,5)
ax2.legend(loc='center left')
plt.title("Turbine in/out pressure with Ambient Temperature (T_des=25°c)")
#plt.title("BH Pressure inlet variation with Ambient Temperature,Florence(T_des=1°c)")
#plt.title("BH Pressure inlet variation with Ambient Temperature,Cheyenne(T_des=1°c)")
#plt.title("BH Pressure inlet variation with Ambient Temperature,Riyadh(T_des=1°c)")

plt.show()
# %%
fig_5_1=plt.plot(annual_T_list,Turbine_eta_OD,color= 'r')
plt.xlabel("Ambient Temperature/ °c")
plt.ylabel("Turbine Efficiency/ % ")
plt.title("Turbine Efficiency variation with Ambient Temperature(T_des=25°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Florence(T_des=24°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Cheyenne(T_des=24°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Riyadh(T_des=24°c)")
plt.ylim(0.5,0.8)
plt.grid()
plt.show()
# %%----- cycle efficiency
fig_3=plt.plot(annual_T_list,cycle_eta,color= 'r')
plt.xlabel("Ambient Temperature/ °c")
plt.ylabel("Cycle Efficiency/ % ")
plt.title("Cycle Efficiency variation with Ambient Temperature (T_des=25°c)")
plt.ylim(0,50)
#plt.title("Cycle Efficiency variation with Ambient Temperature,Florence(T_des=1°c)")
#plt.title("Cycle Efficiency variation with Ambient Temperature,Cheyenne(T_des=1°c)")
#plt.title("Cycle Efficiency variation with Ambient Temperature,Riyadh(T_des=1°c)")
plt.grid()
plt.show()
# %%----- flow coefficient
fig_5=plt.plot(Turbine_m_dot_OD,Turbine_eta_OD,color= 'r')
plt.xlabel("Mass flowrate/ kg/s")
plt.ylabel("Turbine Efficiency/ % ")
plt.title("Turbine Efficiency variation with Mass flowrate,Munich(T_des=1°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Florence(T_des=1°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Cheyenne(T_des=1°c)")
#plt.title("Turbine Efficiency variation with Mass flowrate,Riyadh(T_des=1°c)")
plt.ylim(0.66,0.73)
plt.xlim(22,70)
plt.grid()
plt.show()

# %%----- flow coefficient
max_fi=max(fi_in)
fi_1=list()
for i in fi_in:
    fi_1.append(i/max_fi)
fig_6=plt.plot(fi_in,pressure_ratio,color= 'r')
plt.xlabel("Flow coefficient ")
plt.ylabel("Pressure ratio ")
plt.title("Pressure Ratio variation with Flow Coefficient,Munich(T_des=1°c)")
#plt.title("Pressure Ratio variation with Flow Coefficient,Florence(T_des=1°c)")
#plt.title("Pressure Ratio variation with Flow Coefficient,Cheyenne(T_des=1°c)")
#plt.title("Pressure Ratio variation with Flow Coefficient,Riyadh(T_des=1°c)")
plt.xlim(0.325,1.2)
plt.grid()
plt.show()
# %%----- flow coefficient


fig_7 = plt.plot(Turbine_m_dot_OD,pressure_ratio,color= 'r')
plt.xlabel("Mass flow rate/ kg/s")
plt.ylabel("Pressure Ratio")
plt.title("Pressure Ratio variation with Mass flow rate,Munich(T_des=1°c)")
#plt.title("Pressure Ratio variation with Mass flow rate,Florence(T_des=1°c)")
#plt.title("Pressure Ratio variation with Mass flow rate,Cheyenne(T_des=1°c)")
#plt.title("Pressure Ratio variation with Mass flow rate,Riyadh(T_des=1°c)")
plt.xlim(22,70)
plt.grid()
plt.show()

# %%----- flow coefficient


fig_9 = plt.plot(annual_T_list,pressure_ratio,color= 'r')
plt.ylabel("Pressure ratio")
plt.xlabel("Ambient temperature")
plt.title("Pressure ratio variation with Ambient Temperature,Munich(T_des=1°c)")
#plt.title("Pressure ratio variation with Ambient Temperature,Florence(T_des=1°c)")
#plt.title("Pressure ratio variation with Ambient Temperature,Cheyenne(T_des=1°c)")
#plt.title("Pressure ratio variation with Ambient Temperature,Riyadh(T_des=1°c)")
plt.grid()
plt.show()
