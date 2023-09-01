# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code import constants
import numpy as np
import os


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
point = PlantThermoPoint(["Water"], [1])

t_values = np.linspace(1, 90, 50)   #[°C]
p_values = [0.1, 1, 2, 3, 4]                        #[MPa]
variables = ["Cp", "rho", "k", "mu"]
results = np.zeros((len(t_values), len(p_values), 4))


for i in range(len(t_values)):

    for j in range(len(p_values)):

        point.set_variable("P", p_values[j])
        point.set_variable("T", t_values[i])

        for k in range(len(variables)):
            results[i, j, k] = point.get_variable(variables[k])


# %%------------   EXPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP2", "Task2.2", "4 - Property Comparison", "res", "output"

)

file_path = os.path.join(RES_FOLDER, "properties_UNIFI.xlsx")
units = {

    "Property Units":dict()

}

for k in range(len(variables)):

    data = {

        'T': {"unit": ["°C"], "values": [t_values]},
        'P [MPa]': {"unit": p_values, "values": results[:, :, k].T},

    }

    variable = variables[k]
    write_excel_sheet(excel_path=file_path, sheet_name=variable, data_frame=data, overwrite="hard")
    units["Property Units"].update({variable: {"value": point.get_unit(variable), "unit":None}})

write_excel_sheet(

    excel_path=file_path, sheet_name='Property Units',
    data_frame=units, overwrite="hard",
    write_data=False

)
