from main_code.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.off_design_model.turbine_off_design import TurbineOD
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
# %%--------------------------------------------Read File-------------------------------
Munich=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Munich ERA.csv'
Max_T_amb_Munich =[6.064516129032258, -0.3642857142857144, 6.8903225806451625, 19.220000000000002, 21.25483870967742, 22.803333333333335, 24.819354838709685, 26.348387096774196, 21.07, 15.651612903225804, 7.503333333333334, 4.286666666666665]
Min_T_amb_Munich =[2.9870967741935486, -3.4714285714285715, 0.36129032258064503, 8.02, 11.80967741935484, 13.760000000000002, 14.603225806451615, 15.654838709677417, 10.843333333333332, 7.261290322580646, 2.8866666666666676, 2.4366666666666665]
Max_T_amb_Florence =[11.964516129032264, 7.596428571428572, 12.30645161290323, 21.09333333333333, 23.090322580645164, 27.286666666666676, 31.393548387096775, 31.316129032258065, 28.203333333333333, 22.43870967741935, 15.036666666666664, 10.899999999999999]
Min_T_amb_Florence =[6.425806451612903, 3.0285714285714285, 6.051612903225808, 11.159999999999998, 14.774193548387098, 17.56, 20.690322580645162, 19.99032258064516, 16.63, 13.251612903225805, 9.633333333333333, 4.233333333333333]
Max_T_amb_Riyadh =[21.200000000000006, 25.63214285714286, 31.23870967741936, 32.423333333333325, 38.18709677419355, 42.56333333333334, 43.29677419354839, 42.35483870967742, 41.980000000000004, 33.58387096774193, 25.97666666666667, 23.036666666666665]
Min__T_amb_Riyadh =[8.351612903225806, 13.01785714285714, 17.941935483870974, 20.91, 25.687096774193556, 29.363333333333337, 30.348387096774193, 29.077419354838714, 28.113333333333337, 22.209677419354836, 16.363333333333337, 12.649999999999999]
Max_T_amb_Cheyenne =[3.8258064516129036, 1.8785714285714288, 8.470967741935485, 11.416666666666663, 18.267741935483876, 24.82333333333333, 26.75483870967742, 25.42258064516129, 23.746666666666673, 11.893548387096775, 6.086666666666667, 3.0900000000000003]
Min_T_amb_Cheyenne =[-2.609677419354838, -6.235714285714285, -1.8225806451612903, 0.01666666666666664, 7.474193548387097, 11.753333333333332, 14.151612903225807, 12.919354838709678, 10.150000000000002, 2.335483870967742, -1.3133333333333332, -3.6966666666666668]


# %%--------------------------------------------Power Calculation-------------------------------
def power_off_des(City,t_design):
    Ambient_Temperature = pd.read_csv(City, usecols=['time', 'temperature'])
    Ambient_Temperature['time'] = pd.to_datetime(Ambient_Temperature['time'])
    Monthly_df = \
        Ambient_Temperature[Ambient_Temperature["time"].between
        ('2018-01-01T00:00', '2018-01-31T00:00')]


    annual_T_list = Monthly_df['temperature']
    #Monthly_df["date"] = Monthly_df['time'].dt.date
    #Monthly_df["hour"] = Monthly_df['time'].dt.hour
    #T_amb_avg = Monthly_df.groupby([Monthly_df['time'].dt.month,
                                    #Monthly_df['time'].dt.hour]).temperature.mean()
    powr = tqdm(desc="Calculation", total=len(annual_T_list))
    #annual_T_list=list(T_amb_avg)
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
    heat_sink = list()
    heat_source = list()
    cycle_eta = list()
    delta_T_appr = 5
    Turbine_out_r_points = list()   #turbine out
    eta = 0.7
    T_out = list
    h_out_iso = list()
    h_out_r = list()
    work = list()
    Turbine_Power = list()

    point = PlantThermoPoint(["Carbon Dioxide"], [1])
    for s in annual_T_list:
        T_mean = float(s) + delta_T_appr
        T_Cond_annual.append(T_mean)

    for u in T_Cond_annual:
        tmp_point = point.duplicate()
        tmp_point.set_variable("T", u)
        if u < 30:
            tmp_point.set_variable("Q", 0)
        else:
            tmp_point.set_variable("rho", 700)
        AC_pressure.append(tmp_point.get_variable("P"))
        BH_in_point.append(tmp_point)



    BH_points = PlantThermoPoint(["Carbon Dioxide"], [1])

    dz_well = 1500  # [m] Depth of the reservoir                    #BH output T and P
    T_rock = 125  # [°C] Temperature of the reservoir
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

    TO_point = PlantThermoPoint(["Carbon Dioxide"], [1])
    for iso in range(len(AC_pressure)):
        tmp_point = TO_point.duplicate()
        a = s_in[iso]
        b = AC_pressure[iso]
        tmp_point.set_variable("s", a)
        tmp_point.set_variable("P", b)
        h_out_iso.append(tmp_point.get_variable("h"))
    for j in range(len(AC_pressure)):
        e = h_in[j]
        f = h_out_iso[j]
        x = e - (eta * (e - f))
        h_out_r.append(x)
        work.append(e-x)
    for tor in range(len(AC_pressure)):  # setting real output
        tmp_point = TO_point.duplicate()
        c = h_out_r[tor]
        d = AC_pressure[tor]
        tmp_point.set_variable("h", c)
        tmp_point.set_variable("P", d)
        Turbine_out_r_points.append(tmp_point)

    Turbine_des_m_dot = 10
    for k in range(len(AC_pressure)):
        g = h_in[k] - h_out_r[k]
        work.append(g)
        Turbine_Power.append(float(work[k]) * Turbine_des_m_dot)


    ac_point = PlantThermoPoint(["Carbon Dioxide"], [1])
    for y in T_Cond_annual:
        tmp_point = ac_point.duplicate()
        tmp_point.set_variable("T", y)
        if y < 30:
            tmp_point.set_variable("Q", 0)
        else:
            tmp_point.set_variable("rho", 700)

        h_c_out.append(tmp_point.get_variable("h"))
    for q in range(len(h_out_r)):
        r = h_out_r[q]
        l = h_c_out[q]
        o = h_in[q]
        n = r - l
        u = o - l
        heat_sink.append(n)
        heat_source.append(u)
    for w in range(len(heat_source)):
        x = heat_source[w]
        y = heat_sink[w]
        z = work[w]
        eta = (z / x) * 100
        cycle_eta.append(eta)


    plt.show()

    turb_des_out = PlantThermoPoint(["Carbon Dioxide"], [1])

    t_des= t_design

    turb_des_out.set_variable("T", t_des)

    if t_des < 30:
        turb_des_out.set_variable("Q", 0)
    else:
        turb_des_out.set_variable("rho", 700)


    turb_des_in = PlantThermoPoint(["Carbon Dioxide"], [1])
    dz_well = 1500
    T_rock = 125

    bhe_in = SimplifiedBHE(

        input_thermo_point=turb_des_out,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True)

    bhe_in.update()
    bhe_in.points[-1].copy_state_to(turb_des_in)
    turb_des_in.m_dot = 10

    Turbine_Power_OD = list()
    turbine_d=list()
    turbine = TurbineOD(input_point=turb_des_in, output_point=turb_des_out)

    for j in range(len(BH_in_point)):
        k = BH_in_point[j]
        l = BH_out_points[j]
        turbine.update_input_output_points(l, k)
        turbine.update_off_design_flow_rate()
        Turbine_Power_OD.append(turbine.power)

        powr.update(1)
    powr.close()
    return Turbine_Power_OD

# %%--
def annual_kWh(City,t_design):
    annual_power_production=sum(power_off_des(City,t_design))
    return annual_power_production
annual_kWh(Munich,22)
# %%--
def Sort_data(City):

    Ambient_Temperature = pd.read_csv(City, usecols=['time', 'temperature'])
    Ambient_Temperature['time'] = pd.to_datetime(Ambient_Temperature['time'])
    Monthly_df = \
        Ambient_Temperature[Ambient_Temperature["time"].between
                            ('2018-01-01T00:00', '2018-01-31T00:00')]
    T_amb = Monthly_df['temperature']

    return T_amb

T_amb=Sort_data(Munich)
# %%--

x = power_off_des(Munich,1)
y = power_off_des(Munich,12)

z = power_off_des(Munich,24)



# %%-----------------------------Data Frame-----------------------------------------------------



mean_power = {'T_amb':T_amb,'t_des=1':x,'t_des=12':y,'t_des=24':z}
df_power=pd.DataFrame(mean_power)
df_power.to_csv(r'C:\Users\af_sh\PycharmProjects\power.csv')
print(df_power)
# %%------plotting----

plt.plot(T_amb,x,T_amb,y,T_amb,z)

plt.title('Off-design Power ii January for various T_des, #Munich')
plt.xlabel("T_amb")
plt.ylabel("Power/kW")
plt.legend(('T_des=1','T_des=12','T_des=24'))
plt.show()





