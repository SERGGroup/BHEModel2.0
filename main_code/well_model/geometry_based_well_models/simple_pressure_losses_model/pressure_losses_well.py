from main_code.well_model.simplified_well.simplified_well import SimplifiedWell, SimplifiedBHE, SimplifiedCPG
from main_code.support import PlantThermoPoint
from scipy.constants import g
from abc import ABC
import numpy as np


class PressureLossesWell(SimplifiedWell, ABC):

    def __init__(

            self, input_thermo_point: PlantThermoPoint,
            dz_well, d_inj, d_prod, t_rocks=None, t_surf=None,
            geo_flux=None, k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            heating_section=None, PPI=None, use_rk=True,
            calc_discrete_pressure_losses=False

    ):

        super().__init__(

            input_thermo_point=input_thermo_point, dz_well=dz_well, t_rocks=t_rocks,
            heating_section=heating_section, k_rocks=k_rocks, c_rocks=c_rocks,
            rho_rocks=rho_rocks, t_surf=t_surf, geo_flux=geo_flux, PPI=PPI,
            use_rk=use_rk,

        )

        self.d_inj = d_inj          # unit: m
        self.d_prod = d_prod        # unit: m
        self.calc_discrete_pressure_losses = calc_discrete_pressure_losses

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def additional_final_calculations(self, input_point, depths):

        dh = self.new_point.get_variable("h") - input_point.get_variable("h")
        co_in = self.calculate_C0(input_point, depths[0])
        co_out = self.calculate_C0(self.new_point, depths[1])
        c0 = (co_in + co_out) / 2

        if self.calc_discrete_pressure_losses:

            dp_loss = self.__evaluate_discrete_pressure_losses(input_point, self.new_point, c0=c0)
            self.new_point.set_variable("P", self.new_point.get_variable("P") - dp_loss)
            self.new_point.set_variable("h", input_point.get_variable("h") + dh)

    def evaluate_pressure_losses(self, curr_point: PlantThermoPoint):

        if self.is_upward:
            d = self.d_prod

        else:
            d = self.d_inj

        if self.calc_discrete_pressure_losses or (d is None):

            return 0.

        m_dot = curr_point.m_dot
        mu = curr_point.get_variable("mu")
        rho0 = curr_point.get_variable("rho")

        re = 4 * m_dot / (np.pi * d * mu / 1e6)
        f = self.__calculate_friction_factor(d, re)
        dp = 8 * (f * m_dot ** 2) / (d ** 5 * np.pi ** 2 * rho0)

        return dp / 1e6

    def __evaluate_discrete_pressure_losses(self, input_point:PlantThermoPoint, out_point:PlantThermoPoint, c0):

        m_dot = input_point.m_dot
        mu = (input_point.get_variable("mu") + out_point.get_variable("mu")) / 2
        rho0 = max(input_point.get_variable("rho"), out_point.get_variable("rho"))

        if input_point.get_variable("P") > out_point.get_variable("P"):
            d = self.d_prod

        else:
            d = self.d_inj

        re = 4 * m_dot / (np.pi * d * mu / 1e6)
        f = self.__calculate_friction_factor(d, re)
        l_mod = (np.exp(c0 / 1e6 * g * self.dz_well) - 1) / (g * c0 / 1e6)
        dp = 8 * (f * m_dot ** 2) / (d ** 5 * np.pi ** 2 * rho0) * l_mod

    @staticmethod
    def __calculate_friction_factor(d, re):

        if d is not None:

            rough = 55 * 1e-6       # [m] surface roughness
            rel_rough = rough / d

            A = (2.457 * np.log(1. / ((7. / re) ** 0.9 + 0.27 * rel_rough))) ** 16
            B = (37530. / re) ** 16

            f = 8. * ((8. / re) ** 12 + (1 / (A + B)) ** 1.5) ** (1. / 12.)

        else:

            f = 0.

        return f

    def additional_setup_data(self, data_frame: dict):

        p_losses = "discrete"

        if self.calc_discrete_pressure_losses:

            p_losses = "calculated"

        data_frame["Calculation Options"].update({

            "well class": {"value": "PressureLossesWell", "unit": None},
            "pressure losses": {"value": p_losses, "unit": None}

        })

        return data_frame


class PressureLossesBHE(PressureLossesWell, SimplifiedBHE):

    pass

    def additional_setup_data(self, data_frame: dict):

        data_frame = super().additional_setup_data(data_frame)
        data_frame["Calculation Options"].update({

            "well class": {"value": "PressureLossesBHE", "unit": None},
            "well model": {"value": "BHE", "unit": None}

        })
        return data_frame
