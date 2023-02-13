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
from datetime import datetime, date
# %%--------------------------------------------Read File-------------------------------
Florence=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Florence ERA.csv'
Munich=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Munich ERA.csv'
Cheyenne=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Cheyenne ERA.csv'
Riyadh=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Riyadh ERA.csv'
def Sort_data(City):
    Ambient_Temperature = pd.read_csv(City, usecols=['time', 'temperature'])
    Ambient_Temperature['time'] = pd.to_datetime(Ambient_Temperature['time'])
    Monthly_df = \
        Ambient_Temperature[Ambient_Temperature["time"].between
                            ('2018-01-01T00:00', '2018-12-31T00:00')]

    Monthly_df["date"] = Monthly_df['time'].dt.date
    Monthly_df["hour"] = Monthly_df['time'].dt.hour
    T_amb_avg=Monthly_df.groupby([Monthly_df['time'].dt.month,
                           Monthly_df['time'].dt.hour]).temperature.mean()
    pd.DataFrame(T_amb_avg)

    Mean_T_amb_list = list(T_amb_avg)
    return Mean_T_amb_list  #annual

def Sort_by_month(City):
    monthly_av_list = list()
    b=Sort_data(City)
    for o in range(0, len(b), 24):
        monthly_av_list.append(b[o:o + 24])
    monthly_av_list

    return monthly_av_list
def plot_T_mean(City):
    h = linspace(0, 23, 24)
    for z in range(0, 12):
        month = Sort_by_month(City)[z]
        plt.plot(h, month)
    plt.ylabel("T_amb")
    plt.xlabel("hour")
    plt.legend((

        'Jan', 'Feb', 'Mar',
        'Apr', 'May', 'Jun',
        'Jul', 'Aug', 'Set',
        'Oct', 'Nov', 'Dec'),

        loc='upper right')
    plt.title("hourly T_amb Variation over year 2018")
    plt.show()



# %%-----------------------------------BH input/condenser out point setup-----------------------------------------------------
def Condenser_T(City):
    T_Cond_annual=list()
    delta_T_appr=5
    annual_T_list=Sort_data(City)
    for s in annual_T_list:
        T_mean=float(s)+delta_T_appr
        T_Cond_annual.append(T_mean)
    return T_Cond_annual

# %%-
def Condenser_P(City):
    AC_pressure = list()
    C_point = PlantThermoPoint(["Carbon Dioxide"], [1])
    T_annual = Condenser_T(City)
    for u in T_annual:
        tmp_point = C_point.duplicate()
        tmp_point.set_variable("T", u)
        if u < 30:
            tmp_point.set_variable("Q", 0)
        else:
            tmp_point.set_variable("rho", 700)
        AC_pressure.append(tmp_point.get_variable("P"))
    return AC_pressure

# %%-
def BH_in_points(City):
                          #BH input Pressure list
    BH_in_point = list()
    BH_point = PlantThermoPoint(["Carbon Dioxide"], [1])
    T_annual=Condenser_T(City)
    for ele in T_annual:
        tmp_point = BH_point.duplicate()
        tmp_point.set_variable("T", ele)
        if ele < 30:
            tmp_point.set_variable("Q", 0)
        else:
            tmp_point.set_variable("rho", 700)

        tmp_point.get_variable("P")
        BH_in_point.append(tmp_point)


    return BH_in_point


# %%-
def Turbine_in_points(City):
    h_in=list()
    s_in=list()
    BH_out_points=list()
    dz_well = 1500  # [m] Depth of the reservoir                    #BH output T and P
    T_rock = 125  # [°C] Temperature of the reservoir
    BH_points=PlantThermoPoint(["Carbon Dioxide"], [1])
    BH_in_point=BH_in_points(City)
    for element in BH_in_point:
        bhe_in = SimplifiedBHE(

        input_thermo_point=element,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

        )

        bhe_in.update()
        output_condition = bhe_in.points[-1]
        h_in.append(output_condition.get_variable("h"))
        s_in.append(output_condition.get_variable("s"))
        bhe_in.points[-1].copy_state_to(BH_points)
        BH_out_points.append(BH_points)
    print("h_out_BH",h_in)
    print("s_out_BH", s_in)
    return [BH_out_points,s_in,h_in]

# %%-

def Turbine_out_points(City):
    Turbine_out_r_points=list()
    eta=0.7
    T_out=list
    h_out_iso=list()
    h_out_r=list()

    s_in=Turbine_in_points(City)[1]
    h_in=Turbine_in_points(City)[2]
    P_out=Condenser_P(City)
    TO_point=PlantThermoPoint(["Carbon Dioxide"], [1])
    for iso in range(len(P_out)):
        tmp_point = TO_point.duplicate()
        a=s_in[iso]
        b=P_out[iso]
        tmp_point.set_variable("s",a )
        tmp_point.set_variable("P",b)
        h_out_iso.append(tmp_point.get_variable("h"))
    for i in range(len(P_out)):
        e=h_in[i]
        f=h_out_iso[i]
        x =  e- (eta * (e - f))
        h_out_r.append(x)


    for tor in range(len(P_out)):   #setting real output
        tmp_point = TO_point.duplicate()
        c=h_out_r[tor]
        d=P_out[tor]
        tmp_point.set_variable("h",c)
        tmp_point.set_variable("P",d)
        Turbine_out_r_points.append(tmp_point)

    return [Turbine_out_r_points, h_out_r]

# %%--
def Turbine_Power(City):
    dz_well = 1500  # [m] Depth of the reservoir                    #BH output T and P
    T_rock = 125
    a=Turbine_in_points(City)[2]   #h in and out
    b=Turbine_out_points(City)[1]
    dh=list()
    Turbine_Power=list()
    Turbine_des_m_dot = 10
    for i in range(len(Turbine_out_points(City)[1])):
        g=a[i]-b[i]

        dh.append(g)
        #turbine = TurbineOD(input_point=g, output_point=u)
        #turbine.update_input_output_points(a[i], b[i])
        #turbine.evaluate_design()
        #turbine.update_off_design_flow_rate()
        Turbine_Power.append(float(dh[i])*Turbine_des_m_dot)
    return Turbine_Power
print("Pow",Turbine_Power(Munich))
 # %%-


Power = Turbine_Power(Munich)
monthly_power = list()
for i in range(0, len(Power), 24):
    monthly_power.append(Power[i:i + 24])
print("l", len(monthly_power))
h = linspace(0, 23, 24)
for j in range(0, 12):
    plt.plot(h, monthly_power[j])
plt.legend((

    'Jan', 'Feb', 'Mar',
    'Apr', 'May', 'Jun',
    'Jul', 'Aug', 'Set',
    'Oct', 'Nov', 'Dec'),

    loc='upper right'

)
plt.ylabel("Power/kW")
plt.xlabel("hour")
plt.title("Power Production over the year")
plt.show()














