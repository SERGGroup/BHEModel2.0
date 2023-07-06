# %%------------   IMPORT MODULES                     ----------------------------------------------------------->
from main_code.simplified_well.simplified_well_subclasses import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

# %%------------   BH input condition                   -----------------------------------------------------------> #


class Condenser:
    CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])

    def __init__(self,T_amb:list,dT_appr=7):
        self.T_amb=T_amb
        self.dT_appr=dT_appr
        self.T_c_out=T_c_out
    def condeser_T_out(self):
        T_c_out=self.T_amb+self.dT_appr
        return [T_c_out]

    def condenser_P(self):
        CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
        self.CO2_input.set_variable("Q",0)
        self.CO2_input.set_variable("T",T_c_out)
        P_c = self.CO2_input.get_variable("P")

T1=Condenser([1,2,3],7)
print(T1.T_c_out)
