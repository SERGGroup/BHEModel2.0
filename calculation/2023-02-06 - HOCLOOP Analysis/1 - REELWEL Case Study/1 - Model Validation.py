# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import numpy as np

from main_code.simplified_well.heating_sections.subclasses import REELWELLHeatingSection, REELWELLGeometry
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_well.simplified_well import SimplifiedBHE


# %%------------   IMPORT VALIDATION DATA             -----------------------------------------------------------> #

#   Validation data from:
#
#       "20221107 HOCLOOP - CO2 case well" (in res folder)
#
#   Data extracted from plot using:
#
#       WebPlotDigitizer (https://automeris.io/WebPlotDigitizer)

# TODO

# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

t_in = 30           # [C]
depth = 4500        # [m]
l_overall = 10000   # [m]
r_cur = 500         # [m]

t_rock = 172.5      # [C]
k_rock = 2.44       # [W/(m K)]
c_rock = 1.085      # [kJ/(kg K)]
rho_rock = 2542     # [kg/m^3]

l_horiz = l_overall - (depth + (np.pi / 2 - 1) * r_cur)
hs_geometry = REELWELLGeometry(l_horiz)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #

bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 0.2)
bhe_in.m_dot = 10.

well = SimplifiedBHE(

    bhe_in, dz_well=depth, T_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock

)

well.update()
print(well)

# %%------------   CALCULATIONS                           -----------------------------------------------------------> #

heating_section = REELWELLHeatingSection(well, hs_geometry, hot_in_tubing=False)
well.heating_section = heating_section
well.update()
print(well)
