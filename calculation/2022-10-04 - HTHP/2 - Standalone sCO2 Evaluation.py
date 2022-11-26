# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.power_plants.HTHP.subclasses.CO2_heat_pump_py import StandaloneCO2HeatPump
from main_code.support.other.support_functions import get_np_array
from main_code.support.other.matplolib_stiles import ColorFader
from matplotlib.ticker import MultipleLocator
from main_code import constants
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os

# %%------------   INITIALIZATION AND FIRST TESTS         -----------------------------------------------------------> #

hthp = StandaloneCO2HeatPump(

    BHE_depth=2500,
    T_rock=130,
    P_steam=1,
    T_ambient=15,

)

hthp.thermo_analysis()
print(hthp)

# %%------------   CALCULATIONS                           -----------------------------------------------------------> #

T_SG_perc_list = get_np_array(0.2, 1, 5)
rho_in_list = get_np_array(500, 800, 20)

hthp_list = list()
pbar = tqdm(desc="calculation", total=len(T_SG_perc_list) * len(rho_in_list))

for T_SG_perc in T_SG_perc_list:

    hthp_sublist = list()

    for rho_in in rho_in_list:

        try:

            hthp = StandaloneCO2HeatPump(

                BHE_depth=800,
                T_rock=90,
                P_steam=1,
                T_in_BHE=35,
                T_ambient=15,

            )

            hthp.T_SG_perc = T_SG_perc
            hthp.rho_in = rho_in

            hthp.thermo_analysis()

            pbar.update(1)

        except Exception as error:

            pass

        else:

            hthp_sublist.append(hthp)

    hthp_list.append(hthp_sublist)

pbar.close()

# %%------------   DATA PLOT                              -----------------------------------------------------------> #

# Figure initialization
fig, (ax_1, ax_2) = plt.subplots(1, 2, dpi=150)
fig.set_size_inches(12, 5)

# Math Test Font Definition
plt.rcParams.update({'mathtext.fontset': 'dejavuserif'})

# Font Calculation
T_SG_min = np.min(T_SG_perc_list)
T_SG_max = np.max(T_SG_perc_list)
cf = ColorFader(

    from_value=T_SG_min,
    to_value=T_SG_max

)

for hthp_sublist in hthp_list:

    x_list = list()
    y_list = list()
    y_2_list = list()

    for hthp in hthp_sublist:

        x_list.append(hthp.rho_in)
        y_list.append(hthp.m_dot_ratio)
        y_2_list.append(hthp.m_dot_turb_LP / hthp.m_BHE * 100)

    ax_1.plot(

        x_list, y_list,
        label=r'$T_{{SG\%}}={:.2f}$'.format(hthp.T_SG_perc),
        color=cf.get_color(hthp.T_SG_perc)

    )

    ax_2.plot(

        x_list, y_2_list,
        label=r'$T_{{SG\%}}={:.2f}$'.format(hthp.T_SG_perc),
        color=cf.get_color(hthp.T_SG_perc)

    )

plt.legend(loc='upper right')

y_labels = [

    r'$m_{{ratio}} [-]$',
    r'$m_{{Turb_{{LP}}\%}} [\%]$'

]

i = 0

for ax in [ax_1, ax_2]:

    # Axes Title
    ax.set_xlabel(

        r'$rho_{{in}} [kg/m^3]$',
        fontsize='large', loc='right',
        weight='light'

    )

    ax.set_ylabel(

        y_labels[i],
        fontsize='large', loc='top',
        weight='light'

    )

    i += 1

    # Axis Styling
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax.tick_params(axis=u'both', which=u'minor', length=0)

    # Grid Definition
    ax.set_xlim([500, 800])

    xticks_int = 100
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(xticks_int))
    ax.xaxis.set_minor_locator(MultipleLocator(xticks_int/5))
    ax.grid(which='major', linewidth=0.8)
    ax.grid(which='minor', linewidth=0.2, alpha=0.75)

plt.show()

# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

current_folder = os.path.join(os.path.dirname(constants.RES_FOLDER), "2022-10-04 - HTHP")
output_directory = os.path.join(current_folder, "outputs")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "StandaloneHTHP_behaviour.png"))
