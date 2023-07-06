from main_code.simplified_well.heating_sections.subclasses.REELWELL_heating_section import REELWELLGeometry
from main_code.simplified_well.simplified_well import SimplifiedWell
from main_code.support import PlantThermoPoint


class SimplifiedBHE(SimplifiedWell):

    def evaluate_points(self):

        self.C0[0], self.P_loss[0] = self.update_DP_vertical(self.points[0], is_upward=False)
        self.heating_section.update()
        self.C0[1], self.P_loss[1] = self.update_DP_vertical(self.points[2], is_upward=True)