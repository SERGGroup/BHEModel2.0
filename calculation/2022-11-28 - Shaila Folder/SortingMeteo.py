# %%-----------------------------------Path Setting-------------------------------------------------------
from main_code.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.off_design_model.turbine_off_design import TurbineOD
from pathlib import Path
import pandas as pd
import numpy as np
from numpy import linspace
from pandas import Series,DataFrame
import matplotlib.pyplot as plt
import calendar
# %%--
#Florence
Florence=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Florence ERA.csv'

#Munich
Munich=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Munich ERA.csv'


Ambient_Temperature = pd.read_csv(Munich,usecols=['time','temperature'])
Ambient_Temperature['time']=pd.to_datetime(Ambient_Temperature['time'])
print(Ambient_Temperature.shape)

Monthly_df=Ambient_Temperature[Ambient_Temperature["time"].between('2018-01-01T00:00','2018-12-31T00:00')]
Monthly_df["month"] = Monthly_df['time'].dt.month
Monthly_df["date"] = Monthly_df['time'].dt.date
Monthly_df["hour"] = Monthly_df['time'].dt.hour
print(Monthly_df)
T_amb_avg=(Monthly_df.groupby([Monthly_df['time'].dt.month, Monthly_df['time'].dt.hour]).temperature.mean())
pd.DataFrame(T_amb_avg)
print(T_amb_avg)


# %%-----------------------------------listing-----------------------------------------------------
T_amb_list=list(T_amb_avg)
print(len(T_amb_list))
monthly_av_list=list()
for i in range(0, len(T_amb_list), 24):
    monthly_av_list.append(T_amb_list[i:i + 24])
print(monthly_av_list)


# %%-
January=monthly_av_list[0]
February=monthly_av_list[1]
March=monthly_av_list[2]
April=monthly_av_list[3]
May=monthly_av_list[4]
June=monthly_av_list[5]
July=monthly_av_list[6]
August=monthly_av_list[7]
September=monthly_av_list[8]
October=monthly_av_list[9]
November=monthly_av_list[10]
December=monthly_av_list[11]
m=[]
for i in range(1,13):
    m.append(calendar.month_name[i])

print(m)

# %%-----------------------------------Hourly Temperature -----------------------------------------------------
h=linspace(0,23,24)
for i in range(0,12):
    plt.plot(h, monthly_av_list[i])

plt.title("Florence hourly T_amb Variation over year 2018")
plt.legend((

    'Jan', 'Feb', 'Mar',
    'Apr', 'May', 'Jun',
    'Jul', 'Aug', 'Set',
    'Oct', 'Nov', 'Dec'),

    loc='upper right'

)
plt.ylabel("T_amb")
plt.xlabel("hour")
plt.show()


# %%-----------------------------------BH input point setup-----------------------------------------------------
delta_T_appr=5
T_c=list()
for T_amb_Jan in January:                       #BH input temperature list
    value = T_amb_Jan + delta_T_appr
    T_c.append(value)
print("T_c",T_c)

AC_pressure = list()                            #BH input Pressure list
AC_point = PlantThermoPoint(["Carbon Dioxide"], [1])
BH_in_point=list()

for element in T_c:
    a = float(element)
    tmp_point = AC_point.duplicate()
    tmp_point.set_variable("T", a)
    if a < 30:
        tmp_point.set_variable("Q", 0)
    else:
        tmp_point.set_variable("rho", 700)
    AC_pressure.append(tmp_point.get_variable("P"))
    BH_in_point.append(tmp_point)

print("P_C",AC_pressure)
# %%-----------------------------------BH output/turbine in condition settup-----------------------------------------------------

dz_well = 1500  # [m] Depth of the reservoir                    #BH output T and P
T_rock = 125    # [°C] Temperature of the reservoir

bh_out_point_list = list()
Turbine_in_P=list()
Turbine_in_T=list()
Turbine_in_h=list()
Turbine_in_s=list()
Turbine_in_rho=list()

for bh_in in BH_in_point:

    bhe_in = SimplifiedBHE(

        input_thermo_point=bh_in,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

    )

    bhe_in.update()

    output_condition = bhe_in.points[-1]
    bh_out_point_list.append(output_condition)
    Turbine_in_P.append(output_condition.get_variable("P"))
    Turbine_in_T.append(output_condition.get_variable("T"))
    Turbine_in_rho.append(output_condition.get_variable("rho"))
    Turbine_in_s.append(output_condition.get_variable("s"))
    Turbine_in_h.append(output_condition.get_variable("h"))
print("T_in",Turbine_in_T)




# %%-----------------------------------Set design condition----------------------------------------------------------
Turbine_des_out=list()
Turbine_des_in=list()
Turbine_Power=list()
for i in T_c:
    turb_des_in = PlantThermoPoint(["Carbon Dioxide"], [1])
    turb_des_out = PlantThermoPoint(["Carbon Dioxide"], [1])

    turb_des_out.set_variable("T", i)

    if i < 30:
        turb_des_out.set_variable("Q", 0)
    else:
        turb_des_out.set_variable("rho", 700)

    bhe_in = SimplifiedBHE(

        input_thermo_point=turb_des_out,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

    )

    bhe_in.update()
    bhe_in.points[-1].copy_state_to(turb_des_in)
    turb_des_in.m_dot = 10
    turbine=TurbineOD(input_point=turb_des_in,output_point=turb_des_out)

for i in range(len(BH_in_point)):

    k = BH_in_point[i]
    l = bh_out_list[i]

    turbine.update_input_output_points(l, k)
    turbine.update_off_design_flow_rate()
    Turbine_Power.append(turbine.power)

print(Turbine_Power)



# %%-----------------------------------Set design condition----------------------------------------------------------

T_amb_max=max(T_amb_list)
T_amb_min=min(T_amb_list)
x = 0.5

t_amb_des = (T_amb_max-T_amb_min)*x +T_amb_min
t_des = t_amb_des + delta_T_appr

turb_des_in = PlantThermoPoint(["Carbon Dioxide"], [1])
turb_des_out = PlantThermoPoint(["Carbon Dioxide"], [1])

turb_des_out.set_variable("T", t_des)

if t_des < 30:
    turb_des_out.set_variable("Q", 0)
else:
    turb_des_out.set_variable("rho", 700)

bhe_in = SimplifiedBHE(

    input_thermo_point=turb_des_out,
    dz_well=dz_well, T_rocks=T_rock, use_rk=True

)

bhe_in.update()
bhe_in.points[-1].copy_state_to(turb_des_in)
turb_des_in.m_dot = 10
turbine=TurbineOD(input_point=turb_des_in,output_point=turb_des_out)


# %%-----------------------------------Power Calculation January-----------------------------------------------------
T_c_January=list()

T_turb_out_Jan=list()

for T in January:
    T_c_Jan=T+delta_T_appr
    T_c_January.append(T_c_Jan)
for item in T_c_January:
    f=item * 1.1
    T_turb_out_Jan.append(f)
print(T_turb_out_Jan)
print(T_c_January)



P_c_January=list()
BH_in_Jan_point=list()                      #this is the BH in point list
Jan_point = PlantThermoPoint(["Carbon Dioxide"], [1])
for T_jan in T_c_January:
    T_jan = float(T_jan)
    tmp_point = Jan_point.duplicate()
    tmp_point.set_variable("T", T_jan)
    if T_jan < 30:
        tmp_point.set_variable("Q", 0)
    else:
        tmp_point.set_variable("rho", 700)
    P_c_January.append(tmp_point.get_variable("P"))
    BH_in_Jan_point.append(tmp_point)

#out put BH
dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir

BH_out_Jan_point = list()       #BH output point list
Turbine_in_Jan=list()
Turbine_in_Jan_point = PlantThermoPoint(["Carbon Dioxide"], [1])
for bh_in in BH_in_Jan_point:

    bhe_in = SimplifiedBHE(

        input_thermo_point=bh_in,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

    )

    bhe_in.update()

    output_condition = bhe_in.points[-1]
    BH_out_Jan_point.append(output_condition)
    bhe_in.points[-1].copy_state_to(Turbine_in_Jan_point)
    Turbine_in_Jan_point.m_dot = 10
    Turbine_in_Jan.append(Turbine_in_Jan_point)
Turb_out_point_Jan=list()
Turb_out_Jan=PlantThermoPoint(["Carbon Dioxide"], [1])

for j in range(len(P_c_January)):
    u=P_c_January[j]
    v=T_turb_out_Jan[j]
    tmp_point = Turb_out_Jan.duplicate()
    tmp_point.set_variable("P", u)
    tmp_point.set_variable("T",v)
    Turb_out_point_Jan.append(tmp_point)


# %%--

Turbine_in_Jan_point.m_dot = 10
for o in range(len(P_c_January)):
    k = Turbine_in_Jan[o]
    l = Turb_out_point_Jan[o]
    turbine = TurbineOD(input_point=k, output_point=l)
    turbine.update_input_output_points(k, l)
    turbine.update_off_design_flow_rate()
    Turbine_Power.append(turbine.power)

print(Turbine_Power)


plt.plot(h,Turbine_Power,"b")
plt.show()

Turbine_Power=list()






