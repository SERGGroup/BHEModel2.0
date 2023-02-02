# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.error_band import draw_error_band
from main_code.off_design_model import TurbineOD
import matplotlib.pyplot as plt
from main_code import constants
import numpy as np
import os


# %%------------   INIT TURBINE DESIGN                    -----------------------------------------------------------> #

turbine_in = PlantThermoPoint(["Carbon Dioxide"], [1])
turbine_out = PlantThermoPoint(["Carbon Dioxide"], [1])

#   Design data from:
#
#       Shi, D.; Zhang, L.; Xie, Y.; Zhang, D.
#       "Aerodynamic Design and Off-design Performance Analysis of a Multi-Stage
#       S-CO2 Axial Turbine Based on Solar Power Generation System."
#       doi: https://doi.org/10.3390/app9040714
#       (in res folder)

p_in = 15
t_in = 500
m_dot_in = 170
p_out = 9
n_stage = 3

turbine_in.m_dot = m_dot_in
turbine_in.set_variable("P", p_in)
turbine_in.set_variable("T", t_in)
turbine_out.set_to_expansion_result(p_out, 0.92, turbine_in)

turbine = TurbineOD(turbine_in, turbine_out, n_stages=n_stage, eta_des=0.92)
print(turbine)

# %%------------   INITIALIZE VALIDATION DATA             -----------------------------------------------------------> #

#   Validation data from:
#
#       Shi, D.; Zhang, L.; Xie, Y.; Zhang, D.
#       "Aerodynamic Design and Off-design Performance Analysis of a Multi-Stage
#       S-CO2 Axial Turbine Based on Solar Power Generation System."
#       doi: https://doi.org/10.3390/app9040714
#       (in res folder)
#
#   Data extracted from plot using:
#
#       WebPlotDigitizer (https://automeris.io/WebPlotDigitizer)

POWER_DATA = {

    "power": [

        16.170055452865064, 14.70609981515712, 11.112754158964881,
        10.280961182994455, 9.016635859519410, 6.3881700554528660,
        3.8262476894639565, 1.996303142329019, 1.1977818853974114,
        0.2994454713493511

    ],
    "p_out": [

        7.503748125937031, 7.995502248875562, 8.979010494752625,
        9.494752623688155, 9.986506746626686, 10.970014992503748,
        11.97751124437781, 12.99700149925037, 13.488755622188904,
        13.98050974512743

    ]

}
EFFICIENCY_DATA = {

    "m_dot": [

        93.298969072164950, 109.79381443298970, 124.43298969072166,
        147.73195876288660, 163.81443298969074, 178.04123711340208,
        184.22680412371136, 188.14432989690724, 197.62886597938146,
        201.13402061855672

    ],
    "eta": [

        0.4907127429805617, 0.8224622030237582, 0.8717062634989201,
        0.9015118790496761, 0.9222462203023759, 0.9196544276457884,
        0.9144708423326134, 0.9105831533477322, 0.9041036717062636,
        0.8976241900647949

    ]

}

# %%------------   TURBINE OFF-DESIGN CALCULATION         -----------------------------------------------------------> #

flow_rate_list = list()
p_out_list = list()
power_list = list()
eta_list = list()

for p_out_bar in range(50, p_in*10, 5):

    p_out = p_out_bar / 10  # bar to MPa
    turbine_out.set_to_expansion_result(p_out, 0.90, turbine_in)
    turbine.update_off_design_flow_rate()

    flow_rate_list.append(turbine.input_point.m_dot)
    p_out_list.append(p_out)
    power_list.append(turbine.power/1e3)
    eta_list.append(turbine.eta_iso)

flow_rate_list = np.array(flow_rate_list)
p_out_list = np.array(p_out_list)
power_list = np.array(power_list)
eta_list = np.array(eta_list)

# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #

fig = plt.figure(dpi=150)
fig.set_size_inches(20, 8)
ax_pow = fig.add_subplot(1, 2, 1)
ax_eta = fig.add_subplot(1, 2, 2)

draw_error_band(ax_pow, p_out_list, power_list, 1, err_fixed=True, alpha=0.4, zorder=5)
draw_error_band(ax_eta, flow_rate_list, eta_list, 0.025, alpha=0.4, zorder=5)

ax_pow.plot(POWER_DATA["p_out"], POWER_DATA["power"], "x", zorder=10)
ax_eta.plot(EFFICIENCY_DATA["m_dot"], EFFICIENCY_DATA["eta"], "x", zorder=10)

ax_pow.set_xlim(7, 15)
ax_pow.set_ylim(0, 18)

ax_eta.set_xlim(90, 210)
ax_eta.set_ylim(0.4, 1)

ax_pow.set_xlabel(r'$P_{out}$ [MPa]', fontsize='large', loc='center')
ax_pow.set_ylabel(r'Power [MW]', fontsize='large', loc='center')
ax_pow.set_title("Power Curve", fontsize='xx-large', loc='center')

ax_eta.set_xlabel(r'mass flow [kg/s]', fontsize='large', loc='center')
ax_eta.set_ylabel(r'$\eta_{iso}$ [-]', fontsize='large', loc='center')
ax_eta.set_title("Efficiency Curve", fontsize='xx-large', loc='center')

plt.show()

# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

current_folder = os.path.join(os.path.dirname(constants.RES_FOLDER), "2023-02-02 - Turbine OD validation")
output_directory = os.path.join(current_folder, "res")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "validation.png"))
