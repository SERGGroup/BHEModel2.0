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

class Anual_Power_Variation:
    def __init__(self):
        self.City=City

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
        monthly_av_list = list()
        for i in range(0, len(Mean_T_amb_list), 24):
            monthly_av_list.append(Mean_T_amb_list[i:i + 24])
        monthly_av_list
                                    #list maming

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
        return monthly_av_list