from main_code.power_plants.subclasses.HTHP.CO2_heat_pump_py import StandaloneCO2HeatPump
from main_code.support.other.support_functions import get_np_array
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# %%-------------------------------------     INITIALIZATION AND FIRST TEST     -------------------------------------> #

hthp = StandaloneCO2HeatPump(

    BHE_depth=3000,
    T_rock=130,
    P_steam=1,
    T_ambient=15,

)

hthp.thermo_analysis()
print(hthp)

# %%-------------------------------------    OPTIMIZATION PARAMETER ANALYSIS    -------------------------------------> #
# %% 1 - Calculations

T_SG_perc_list = get_np_array(0.2, 1, 5)
rho_in_list = get_np_array(550, 800, 20)

hthp_list = list()

for T_SG_perc in T_SG_perc_list:

    hthp_sublist = list()

    for rho_in in rho_in_list:

        try:

            hthp = StandaloneCO2HeatPump(

                BHE_depth=3000,
                T_rock=130,
                P_steam=1,
                T_in_BHE=30,
                T_ambient=15,

            )

            hthp.T_SG_perc = T_SG_perc
            hthp.rho_in = rho_in

            hthp.thermo_analysis()

            print("{:.2f} - {:.2f}".format(T_SG_perc, rho_in))

        except Exception as error:

            pass

        else:

            hthp_sublist.append(hthp)

    hthp_list.append(hthp_sublist)

print("DONE!")

# %% 2 - Data Plot

fig, (ax_1) = plt.subplots(1, 1)
cmap = cm.get_cmap("Blues")

T_SG_min = np.min(T_SG_perc_list)
T_SG_max = np.max(T_SG_perc_list)

for hthp_sublist in hthp_list:

    x_list = list()
    y_list = list()

    for hthp in hthp_sublist:

        x_list.append(hthp.rho_in)
        y_list.append(hthp.points[0].get_variable("P"))

    ax_1.plot(

        x_list, y_list,
        label="T_SG_perc={:.2f}".format(hthp.T_SG_perc),
        color=cmap(hthp.T_SG_perc)

    )

ax_1.set_xlabel("rho_in [kg/m3]")
ax_1.set_ylabel("P_out_BHE [C]")
plt.legend(loc='upper right')
plt.show()
