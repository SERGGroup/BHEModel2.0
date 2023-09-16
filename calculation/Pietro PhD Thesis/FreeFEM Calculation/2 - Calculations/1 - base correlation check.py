# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code.constants import CALCULATION_FOLDER, os
from collections.abc import Iterable
import matplotlib.pyplot as plt
import numpy as np


# %%------------   DEFINE WORKING FOLDERS                 -----------------------------------------------------------> #
calculation_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "calculation folder"

)

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "results"

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
a = 4
DT = 40.
n_mesh = 4
m_mesh = 20
time_steps = 100
n_save = time_steps
results = dict()

mo = MeshOptions(

    n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
    mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
    retrieve_mesh=False, add_visualization_mesh=True

)

pdo = ProblemDefinitionOptions(

    grad_rock=0, DT=DT,
    time_range=[1e-3, 1000], time_steps=time_steps,
    n_plot=0, n_save=n_save

)

tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=False)
ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)

print("CALCULATION STARTED ...")
gross_result = ffa.calculate()
print("DONE!")

time_list = list()
value_list = list()

for line in gross_result:

    line_split = line.split(";")
    time_list.append(float(line_split[0].strip()))
    value_list.append(float(line_split[1].strip()) / (DT * 2 * np.pi))

results.update({

    "time": time_list,
    "value": value_list

})


# %%------------   EVALUATE CORRELATION                   -----------------------------------------------------------> #
def calculate_q(td):

    if issubclass(type(td), Iterable):

        result = list()
        for td_item in td:
            result.append(calculate_q(td_item))

        return result

    else:

        if td < 2.8:

            f = ((np.pi * td) ** -0.5 + 0.5 - (td / np.pi) ** 0.5 / 4 + td / 8)

        else:

            gamma = 0.57722
            theta = (np.log(4 * td) - 2 * gamma)
            f = 2 * (1 / theta - gamma / theta ** 2)

        return f

time_array = np.array(results["time"])
f_array = calculate_q(time_array)


# %%------------   SAVE RESULTS                           -----------------------------------------------------------> #
file_path = os.path.join(result_folder, "1 - base comparison.xlsx")

dataframe = {

    'Time': {"unit": ["-"], "values": [results["time"]]},
    'FreeFEM': {"unit": ["-"], "values": [results["value"]]},
    'Correlation': {"unit": ["-"], "values": [f_array]},

}
write_excel_sheet(excel_path=file_path, sheet_name="results", data_frame=dataframe, overwrite="hard")


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
plt.plot(results["time"], results["value"], label="FreeFEM")
plt.plot(time_array, f_array, "--", label="correlation", color="Black")

plt.xlabel("$t_d$ [-]")
plt.ylabel("$f$ [-]")
plt.xscale("log")
plt.yscale("log")
plt.title("Zhang's Correlation Comparison")
plt.show()
