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
    for i in range(0, len(b), 24):
        monthly_av_list.append(b[i:i + 24])
    monthly_av_list
                                #list maming

    Jan=monthly_av_list[0]
    Feb=monthly_av_list[1]
    Mar=monthly_av_list[2]
    Apr=monthly_av_list[3]
    May=monthly_av_list[4]
    Jun=monthly_av_list[5]
    Jul=monthly_av_list[6]
    Aug=monthly_av_list[7]
    Set=monthly_av_list[8]
    Oct=monthly_av_list[9]
    Nov=monthly_av_list[10]
    Dec=monthly_av_list[11]
    return monthly_av_list #Jan,Feb, Mar, Apr, May, Jun , Jul, Aug, Set, Oct, Nov, Dec
def plot_T_mean(City):
    h = linspace(0, 23, 24)
    for i in range(0, 12):
        month = Sort_by_month(City)[i]
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
Sort_data(Munich)
plot_T_mean(Munich)

# %%-----------------------------------BH input/condenser out point setup-----------------------------------------------------
def Condenser_T(City):
    T_Cond_annual=list()
    delta_T_appr=5
    annual_T_list=Sort_data(City)
    for i in annual_T_list:
        T_mean=float(i)+delta_T_appr
        T_Cond_annual.append(T_mean)
    return T_Cond_annual
Condenser_T(Munich)

# %%-

def BH_in_points(City):
    AC_pressure = list()                            #BH input Pressure list
    BH_in_point = list()
    BH_point = PlantThermoPoint(["Carbon Dioxide"], [1])
    T_annual=Condenser_T(City)
    for i in T_annual:
        tmp_point = BH_point.duplicate()
        tmp_point.set_variable("T", i)
        if i < 30:
            tmp_point.set_variable("Q", 0)
        else:
            tmp_point.set_variable("rho", 700)
        AC_pressure.append(tmp_point.get_variable("P"))

        BH_in_point.append(tmp_point)

    return BH_in_point
def Turbine_in(City):