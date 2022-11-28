# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE

# %%------------   INIT INPUT CONDITION                   -----------------------------------------------------------> #
T_amb = 15  # [°C]
dT_appr = 7  # [°C]

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("Q", 0.)

p_in = CO2_input.get_variable("P")
print(CO2_input.get_unit("D"))


# %%------------   INIT BHE WELL                          -----------------------------------------------------------> #
dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir

bhe_in = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=dz_well, T_rocks=T_rock, use_rk=True

)

bhe_in.update()
print(bhe_in)

output_condition = bhe_in.points[-1]
p_out = output_condition.get_variable("P")
t_out = output_condition.get_variable("T")
h_out = output_condition.get_variable("h")
s_out = output_condition.get_variable("s")

unit_p = output_condition.get_unit("P")

# %%------------   TURBINE EXPANSION                      -----------------------------------------------------------> #
