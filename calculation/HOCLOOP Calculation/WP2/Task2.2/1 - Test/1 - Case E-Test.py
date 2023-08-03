# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import REELWEELBHE, REELWELLGeometry
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
import matplotlib.pyplot as plt
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-a.pdf" (in "0 - Resources\Case Studies" Folder)
#

t_in = 14.72        # [C]
depth = 1828.8      # [m]
mass_flow = 8.8     # [kg/s]

t_surf = 21.111     # [C]
t_grad = 0.01513    # [C/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_surf + t_grad * depth
hs_geometry = REELWELLGeometry(

    depth,
    tub_id=0.05,
    tub_od=0.08,
    cas_id=0.1617,
    cas_od=0.1778,
    k_insulation=0.1,
    hot_in_tubing=True,
    neglect_internal_heat_transfer=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 1)

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    t_surf=t_surf, rw_geometry=hs_geometry, max_iteration=10

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5
step = 50
time = 7

profile_positions = np.linspace(0, depth, int(depth / step + 1))[1:]

well.heating_section.time = time / 365
well.update()


# %%------------   PLOT                                   -----------------------------------------------------------> #
fig, ax = plt.subplots()
profile_pos_ovr = np.concatenate((profile_positions, np.flip(profile_positions)), axis=0)

for profile in well.get_profiles(profile_positions):

    profile_t_ovr = np.concatenate((profile[0][0, :], np.flip(profile[0][1, :])), axis=0)
    profile_p_ovr = np.concatenate((profile[1][0, :], np.flip(profile[1][1, :])), axis=0)
    ax.plot(profile_pos_ovr, profile_t_ovr)

#plt.xscale("log")
plt.show()

print(well.points[-1].get_variable("T"))