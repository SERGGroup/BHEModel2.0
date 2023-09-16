# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt


# %%------------   DEFINE WORKING FOLDERS                 -----------------------------------------------------------> #
calculation_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "calculation folder"

)


# %%------------   INITIALIZATION                         -----------------------------------------------------------> #
a = 4
n_mesh = 2
m_mesh = 5
mo = MeshOptions(

    n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
    mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
    retrieve_mesh=False, add_visualization_mesh=True

)

pdo = ProblemDefinitionOptions(

    grad_rock=0.01, DT=40., pe=1., vx=1, vy=1,
    time_range=[1e-3, 10], time_steps=100,
    n_plot=20, n_save=300

)

tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=False)
ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)
ffa.options.write_edp_file()
#gross_result = ffa.calculate()


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
time_list = list()
value_list = list()

for line in gross_result:

    line_split = line.split(";")
    time_list.append(float(line_split[0].strip()))
    value_list.append(float(line_split[1].strip()))

plt.plot(time_list, value_list, label="{}x - {}%".format(n_mesh, round(100 / m_mesh)))


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
plt.xlabel("$t_d$ [-]")
plt.ylabel("$dT/dr$ [°C/m]")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
