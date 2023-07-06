from main_code.simplified_well.simplified_well import SimplifiedWell
from main_code.support import PlantThermoPoint


class SimplifiedBHE(SimplifiedWell):

    def evaluate_points(self):

        self.C0[0], self.P_loss[0] = self.update_DP_vertical(self.points[0], is_upward=False)
        self.heating_section.update()
        self.C0[1], self.P_loss[1] = self.update_DP_vertical(self.points[2], is_upward=True)


class SimplifiedCPG(SimplifiedWell):

    def __init__(

        self, input_thermo_point: PlantThermoPoint, dz_well, T_rocks,
        heating_section=None, k_rocks=0.2, geo_flux=0.1, q_up_down=0.0,
        PPI=None, use_rk=False, d_inj=None, d_prod=None,
        discretize_p_losses=False, p_res=None

    ):

        super().__init__(

            input_thermo_point, dz_well, T_rocks,
            heating_section=heating_section, k_rocks=k_rocks,
            geo_flux=geo_flux, q_up_down=q_up_down,
            PPI=PPI, use_rk=use_rk, d_inj=d_inj, d_prod=d_prod,
            discretize_p_losses=discretize_p_losses

        )

        if p_res is not None:

            self.p_res = p_res

        else:

            self.p_res = 1000 * 9.81 * dz_well / 1e6  # water hydrostatic pressure (in MPa)

    def evaluate_points(self):

        # iterate point 0 in order to get fixed point 1 (production well bottom hole) pressure.
        # pressure = reservoir pressure - pressure losses in the reservoir

        counter = 0
        surface_temperature = self.points[0].get_variable("T")

        while True:

            self.C0[0], self.P_loss[0] = self.update_DP_vertical(self.points[0], is_upward=False)
            self.heating_section.update()
            dp = self.points[2].get_variable("P") - self.p_res

            if abs(dp / self.p_res) < 1e-3 or counter > 10:

                break

            else:

                counter += 1
                self.points[0].set_variable("P", self.points[0].get_variable("P") - dp)
                self.points[0].set_variable("T", surface_temperature)

        self.C0[1], self.P_loss[1] = self.update_DP_vertical(self.points[2], is_upward=True)
