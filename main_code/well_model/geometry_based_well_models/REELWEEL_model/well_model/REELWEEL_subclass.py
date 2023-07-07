from main_code.well_model.geometry_based_well_models.REELWEEL_model.heating_sections import REELWELLGeometry
from main_code.well_model.simplified_well.heating_sections import EmptyHeatingSection
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support import PlantThermoPoint
from scipy.constants import g


class REELWEELBHE(SimplifiedBHE):

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well,
            t_rocks=None, t_surf=None, geo_flux=None,
            k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            heating_section=None, PPI=None,
            use_rk=True, rw_geometry=None

    ):

        if heating_section is None:

            heating_section = EmptyHeatingSection(self)

        super().__init__(

            input_thermo_point, dz_well, t_rocks=t_rocks,
            heating_section=heating_section, k_rocks = k_rocks, c_rocks = c_rocks, rho_rocks = rho_rocks,
            t_surf=t_surf, PPI = PPI, use_rk = use_rk, geo_flux=geo_flux

        )

        if rw_geometry is not None:

            self.rw_geometry = rw_geometry

        else:

            self.rw_geometry = REELWELLGeometry(dz_well)

        self.rw_geometry.parent_class = self.heating_section

    def evaluate_pressure_losses(self, curr_point: PlantThermoPoint):

        dp = self.rw_geometry.dp_dl(curr_point, is_annulus = self.rw_geometry.hot_in_tubing ^ self.is_upward)
        return dp

    def dh_dl_stream(self, curr_point: PlantThermoPoint, depth=0.):

        is_annulus = self.rw_geometry.hot_in_tubing ^ self.is_upward
        dh = self.rw_geometry.dh_dl(curr_point, is_annulus=is_annulus, depth=depth)
        (g + dh * 1e3) / g * 1e3

        rho = curr_point.get_variable("rho")
        return (g + dh * 1e3) / (g * rho) * 1e3

    def additional_setup_data(self, data_frame: dict):

        data_frame["Calculation Options"].update({

            "well class": {"value": "REELWEEL_BHE", "unit": None}

        })

        data_frame = self.rw_geometry.additional_setup_data(data_frame)
        return data_frame
