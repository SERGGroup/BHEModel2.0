from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.power_plants.off_design_model.turbine_off_design import TurbineOD
from pathlib import Path
import pandas as pd
import numpy as np
from numpy import linspace
from pandas import Series,DataFrame
import matplotlib.pyplot as plt
import calendar
from tqdm import tqdm
from datetime import datetime, date

t_des = 20
tmp_point = BH_point.duplicate()
tmp_point.set_variable("T", ele)

if ele < 30:

    tmp_point.set_variable("Q", 0)

else:

    tmp_point.set_variable("rho", 700)

tmp_point.get_variable("P")
BH_in_point.append(tmp_point)

turbine = TurbineOD(input_point=input_point, output_point=output_point)

temperature_list = [15, 16, 17, 20, ...] #this should czame from the ambient data

for T in range(len(temperature_list)):

    # evaluate input and outpu of the turbine

    turbine.update_input_output_points(new_input_point, new_output_point)
    turbine.update_off_design_flow_rate()

    Turbine_Power.append(turbine.power)
    pbar.update(1)