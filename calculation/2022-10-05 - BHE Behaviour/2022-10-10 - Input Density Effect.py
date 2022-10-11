# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.support_functions import get_np_array
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE
import matplotlib.pyplot as plt
from tqdm import tqdm


# %%-------------------------------------   INITIALIZATION                      -------------------------------------> #

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", 35)
CO2_input.set_variable("rho", 800)

bhe_CO2 = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=1200, T_rocks=150

)

n_points = 100
rho_in_points = get_np_array(800, 1200, n_points)    # in m

# %%-------------------------------------   CO2 CALCULATION                     -------------------------------------> #

pbar_rho = tqdm(desc="evaluating rho analysis:", total=n_points)

x_list = list()
y_list = list()

for rho in rho_in_points:

    CO2_input.set_variable("rho", rho)
    bhe_CO2.update()

    if 0 <= bhe_CO2.eta_II <= 1:

        x_list.append(rho)
        y_list.append(bhe_CO2.eta_II)

    pbar_rho.update(1)

pbar_rho.close()

# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, (ax) = plt.subplots(1, 1, dpi=200)
fig.set_size_inches(10, 5)

ax.plot(

    x_list, y_list,
    linewidth=2

)

plt.show()