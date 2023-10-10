from REFPROPConnector import ThermodynamicPoint
#%%
point1=ThermodynamicPoint(["CarbonDioxide"],[0.])
point1.set_variable("s",1.82)
point1.set_variable("P",13.16)
h_iso=point1.get_variable("h")
print(h_iso)

bhe_in2 = SimplifiedBHE(

    input_thermo_point=CO2_input2,
    dz_well=dz_well, T_rocks=T_rock, use_rk=True

)

bhe_in2.update()
print(bhe_in2)


output_condition = bhe_in2.points[-1]
p_out2 = output_condition.get_variable("P")
t_out2 = output_condition.get_variable("T")
h_out2 = output_condition.get_variable("h")
s_out2 = output_condition.get_variable("s")
unit_p = output_condition.get_unit("P")





T_amb2 = 22 # [°C]
dT_appr = 7  # [°C]

CO2_input2 = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input2.set_variable("T", T_amb2 + dT_appr)
CO2_input2.set_variable("Q",0.)

p_in2 = CO2_input2.get_variable("P")
print(CO2_input2.get_unit("D"))

h_t_in=Turbine_input_condition.get_variable("h")
s_t_in=Turbine_input_condition.get_variable("s")


print(h_out_iso)