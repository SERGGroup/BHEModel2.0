# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.other.excel_exporter import export_profiles_to_excel
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code import constants
import matplotlib.pyplot as plt
import numpy as np
import os


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

t_in = 45           # [C]
depth = 3000        # [m]
l_overall = 6500    # [m]
mass_flow = 15.80149 # [kg/s]

t_surf = 11         # [C]
t_grad = 0.0325     # [c/m]
k_rock = 2.423      # [W/(m K)]
c_rock = 0.90267    # [kJ/(kg K)]
rho_rock = 2600     # [kg/m^3]

t_rock = t_surf + t_grad * depth
l_horiz = l_overall - depth
hs_geometry = REELWELLGeometry(

    l_horiz,
    tub_id=0.1,
    tub_od=0.13,
    cas_id=0.1617,
    cas_od=0.1778,
    k_insulation=0.1,
    hot_in_tubing=True,
    neglect_internal_heat_transfer=False,
    max_back_time=2,
    alpha_old=0.5

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 1)
bhe_in.m_dot = mass_flow

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    rw_geometry=hs_geometry, max_iteration=50

)

heating_section = REELWELLHeatingSection(

    well, hs_geometry,
    integrate_temperature=True,

)
well.heating_section = heating_section


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
n_points = 30
shape = 2.5
step = 5
time = 7

profile_positions = np.linspace(0, depth, int(depth / step + 1))[1:]
hs_profile_positions = np.linspace(0, l_horiz, int(l_horiz / step + 1))[1:]

well.heating_section.time = time / 365
well.update()


# %%------------   PLOT                                   -----------------------------------------------------------> #
fig, ax = plt.subplots()
profile_pos_ovr = np.concatenate((profile_positions, np.flip(profile_positions)), axis=0)

for profile in well.heating_section.get_profiles(profile_positions)[6:7]:

    profile_t_ovr = np.concatenate((profile[0][0, :], np.flip(profile[0][1, :])), axis=0)
    profile_p_ovr = np.concatenate((profile[1][0, :], np.flip(profile[1][1, :])), axis=0)
    ax.plot(profile_pos_ovr, profile_t_ovr)

#plt.xscale("log")
plt.show()

print(well.points[-1].get_variable("T"))


# %%------------   CHECK                                  -----------------------------------------------------------> #
profile = well.heating_section.get_profiles(hs_profile_positions)[-1]


# %%------------   CHECK                                  -----------------------------------------------------------> #
fig, ax = plt.subplots()
dx = hs_profile_positions[1] - hs_profile_positions[0]
r_res = list()
for i in range(len(profile[0][0, :]) - 1):

    dt_in = profile[0][0, i] - profile[0][1, i]
    dt_out = profile[0][0, i+1] - profile[0][1, i+1]
    dtlm = (dt_in - dt_out)/np.log(dt_in / dt_out)
    dh_tub = profile[3][1, i] - profile[3][1, i+1]
    q = well.m_dot * dh_tub
    r_res.append(dtlm / (q / dx * 1000))

ax.plot(hs_profile_positions[:-2], np.array(r_res[:-1]) / hs_geometry.r_ins)

#plt.xscale("log")
plt.show()
print(well.m_dot)
print(hs_geometry.r_ins/np.array(r_res[:-1]))