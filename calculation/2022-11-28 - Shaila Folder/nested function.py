# %%
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE
from tabulate import tabulate
import matplotlib

# %%
delta_pump_p=0.7
CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
class Ambient_T:

  def __iter__(self):
    self.a = 1
    return self

  def __next__(self):

    if self.a <= 42:
      x = self.a
      self.a += 1


      CO2_input.set_variable("T", x)
      if x < 30:


        CO2_input.set_variable("Q", 0.)
        p_c = CO2_input.get_variable("P")  # condenser pressure or turbine outlet pressure

      else:

        CO2_input.set_variable("rho", 700)
        p_c = CO2_input.get_variable("P")

      p_c = CO2_input.get_variable("P")
      print(p_c)
      #CO2_input.set_variable("T", x)
      p_in_BH = p_c + delta_pump_p
      CO2_input.set_variable("P", p_in_BH)
      h_in_BH = CO2_input.get_variable("h")

      return CO2_input
    else:
      raise StopIteration

Temperature_Range = Ambient_T()
T_amb = iter(Temperature_Range)

for x in Temperature_Range:
  print(CO2_input)


  # %%------------   INIT BHE WELL                          -----------------------------------------------------------> #








