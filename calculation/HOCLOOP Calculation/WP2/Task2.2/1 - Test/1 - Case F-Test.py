# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

    REELWELLHeatingSection,
    REELWELLGeometry,
    REELWEELBHE,

)
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
import scipy.constants
import numpy as np


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
#   Validation data from:
#       "slides-case-2.2-b.pdf" (in "0 - Resources\Case Studies" Folder)
#

t_in = 45           # [C]
depth = 3000        # [m]
l_overall = 6500    # [m]
mass_flow = 8.80149 # [kg/s]

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
    max_back_time=4,
    alpha_old=0.5,
    neglect_internal_heat_transfer=False,
    ignore_tubing_pressure_losses=False,
    ignore_annulus_pressure_losses=False

)


# %%------------   INITIALIZE WELL                        -----------------------------------------------------------> #
bhe_in = PlantThermoPoint(["Water"], [1])
bhe_in.set_variable("T", t_in)
bhe_in.set_variable("P", 1.2)
bhe_in.m_dot = mass_flow

well = REELWEELBHE(

    bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=t_surf,
    k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
    rw_geometry=hs_geometry, max_iteration=30

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


# %%------------   PLOT BASIC PROFILES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
profile_pos_ovr = np.concatenate((profile_positions, np.flip(profile_positions)), axis=0)

for profile in well.heating_section.get_profiles(profile_positions)[7:]:

    profile_t_ovr = np.concatenate((profile[0][0, :], np.flip(profile[0][1, :])), axis=0)
    profile_p_ovr = np.concatenate((profile[1][0, :], np.flip(profile[1][1, :])), axis=0)
    ax.plot(profile_pos_ovr, profile_t_ovr)

#plt.xscale("log")
plt.show()

print(well.points[-1].get_variable("T"))


# %%------------   CHECK HEAT TRANSFER                    -----------------------------------------------------------> #
profile = well.heating_section.get_profiles(position_list=hs_profile_positions, get_index=-1)

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


# %%------------   CHECK INTERNAL SPEAD                   -----------------------------------------------------------> #
well_profile = well.get_profiles(position_list=profile_positions, get_index=-1)
hs_profile = well.heating_section.get_profiles(position_list=hs_profile_positions, get_index=-1)

fluid_speed = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
reynolds = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
viscosity = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
friction_factor = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
dp_dl = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
p_profile = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
dp_calc = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
dp_grav = np.zeros((2, len(profile_positions) + len(hs_profile_positions)))
__tmp_point = bhe_in.duplicate()

well_positions = np.concatenate((

    profile_positions,
    hs_profile_positions + profile_positions.max()

))

p_real = np.concatenate((

    profile_positions,
    hs_profile_positions + profile_positions.max()

))

i_0 = 0
l_0 = 0.
profiles = [well_profile, hs_profile]

for profile in profiles:

    profile_len = len(profile[0][0, :])

    for k in range(profile_len):

        i = i_0 + k

        for j in range(2):

            is_annulus = (j == 0)

            __tmp_point.set_variable("P", profile[1][j, k])
            __tmp_point.set_variable("rho", profile[2][j, k])

            fluid_speed[j, i] = hs_geometry.fluid_speed(is_annulus=is_annulus, point=__tmp_point)
            friction_factor[j, i] = hs_geometry.f(is_annulus=is_annulus, point=__tmp_point)
            reynolds[j, i] = hs_geometry.re(is_annulus=is_annulus, point=__tmp_point)
            dp_dl[j, i] = - hs_geometry.dp_dl(is_annulus=is_annulus, point=__tmp_point) * 1e6

            viscosity[j, i] = __tmp_point.get_variable("mu")
            p_profile[j, i] = __tmp_point.get_variable("P")
            dp_calc[j, i] = dp_dl[j, i] * step / 1e6 + dp_calc[j, i-1]

            if i_0 == 0:
                dp_grav[j, i] += (__tmp_point.get_variable("rho") * scipy.constants.g * step) / 1e6 + dp_grav[j, i-1]
            else:
                dp_grav[j, i] += 0. + dp_grav[j, i - 1]

        __tmp_point.set_variable("T", profile[0][1, k])
        __tmp_point.set_variable("rho", profile[2][1, k])

    i_0 += profile_len


# %%------------   PLOT INTERNAL SPEAD                    -----------------------------------------------------------> #
data_to_plot = {

    "Fluid Speed [m/s]": fluid_speed,
    "Friction Factor [-]": friction_factor,
    "Reynolds Number [-]": reynolds,
    "Pressure Losses [Pa/m]": dp_dl,
    "Viscosity [uPa/s]": viscosity

}

for key in data_to_plot.keys():

    fig, ax = plt.subplots()
    ax.plot(well_positions, data_to_plot[key][0, :], label="annulus")
    ax.plot(well_positions, data_to_plot[key][1, :], label="tubing")
    ax.set_xlabel("Position [m]")
    ax.set_ylabel(key)
    ax.legend()

    plt.tight_layout(pad=2)
    plt.show()


# %%------------   PLOT PRESSURE LOSSES                    -----------------------------------------------------------> #
fig, ax = plt.subplots()
ax.plot(well_positions, dp_grav[0, :] - dp_calc[0, :], label="annulus", color="tab:blue")
ax.plot(well_positions, dp_grav[1, :] - dp_calc[1, :], label="tubing", color="tab:orange")

ax.plot(well_positions, p_profile[0, :] - p_profile[0, 0], '--', label="annulus", color="tab:blue")
ax.plot(well_positions, p_profile[1, :] - p_profile[1, 0], '--', label="tubing", color="tab:orange")

ax.set_xlabel("Position [m]")
ax.set_ylabel("Pressure [MPa]")
ax.legend()

plt.tight_layout(pad=2)
plt.show()