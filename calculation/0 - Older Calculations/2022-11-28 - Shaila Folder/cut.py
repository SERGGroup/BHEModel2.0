CO2_input_l = PlantThermoPoint(["CarbonDioxide"], [1])

CO2_input_l.set_variable("T", T_amb_l + dT_appr)
CO2_input_l.set_variable("Q", 0.)
p_in_l = CO2_input_l.get_variable("P")

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("Q", 0.)
p_in = CO2_input.get_variable("P")

CO2_input_h = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input_h.set_variable("T", T_amb_h + dT_appr)
CO2_input_h.set_variable("Q", 0.)
p_in_h = CO2_input_h.get_variable("P")