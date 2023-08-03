from main_code.well_model.geometry_based_well_models.REELWEEL_model.geometry import REELWELLGeometry, ITERATION_TOLLS
from main_code.well_model.simplified_well.heating_sections import EmptyHeatingSection
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support import PlantThermoPoint
from scipy.constants import g
from copy import deepcopy
import numpy as np


class REELWEELBHE(SimplifiedBHE):

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well,
            t_rocks=None, t_surf=None, geo_flux=None,
            k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            heating_section=None, PPI=None,
            use_rk=True, rw_geometry=None,
            max_iteration=10

    ):

        if heating_section is None:

            heating_section = EmptyHeatingSection(self)

        super().__init__(

            input_thermo_point, dz_well, t_rocks=t_rocks,
            heating_section=heating_section, k_rocks=k_rocks, c_rocks=c_rocks, rho_rocks=rho_rocks,
            t_surf=t_surf, PPI=PPI, use_rk=use_rk, geo_flux=geo_flux

        )

        if rw_geometry is not None:

            self.rw_geometry = rw_geometry

        else:

            self.rw_geometry = REELWELLGeometry(dz_well)

        self.__old_profiles = None
        self.__max_iteration = max_iteration

        self.__tmp_point_annulus = input_thermo_point.duplicate()
        self.__tmp_point_tubing = input_thermo_point.duplicate()

        self.rw_geometry.parent_class = self.heating_section

    def evaluate_points(self):

        self.__old_profiles = list()
        self.integrators_profiler = list()
        super().evaluate_points()

        if not self.neglect_internal_heat_transfer:

            for i in range(self.__max_iteration):

                self.__old_profiles.append(deepcopy(self.integrators_profiler))
                self.integrators_profiler = list()
                super().evaluate_points()

                if self.convergence_reached:

                    break

        self.__old_profiles.append(deepcopy(self.integrators_profiler))

    def evaluate_pressure_losses(self, curr_point: PlantThermoPoint):

        dp = self.rw_geometry.dp_dl(curr_point, is_annulus = self.rw_geometry.hot_in_tubing ^ self.is_upward)
        return dp

    def dh_dl_stream(self, curr_point: PlantThermoPoint, depth=0.):

        is_annulus = self.rw_geometry.hot_in_tubing ^ self.is_upward
        dh_ext = self.rw_geometry.dh_dl(curr_point, is_annulus=is_annulus, depth=depth)
        dh_int = self.get_old_dh_int(is_annulus=is_annulus, depth=depth)
        dh = dh_ext + dh_int

        if self.is_upward:
            return g / 1e3 - dh

        else:
            return g / 1e3 + dh

    def additional_setup_data(self, data_frame: dict):

        data_frame["Calculation Options"].update({

            "well class": {"value": "REELWEEL_BHE", "unit": None}

        })

        data_frame = self.rw_geometry.additional_setup_data(data_frame, is_well=True)
        return data_frame

    @property
    def neglect_internal_heat_transfer(self):

        return self.rw_geometry.neglect_internal_heat_transfer

    @property
    def convergence_reached(self):

        if len(self.__old_profiles) > 0:

            check_points = np.linspace(start=0, stop=1, num=50) * self.dz_well
            curr_profiles = self.get_iteration_points(check_points)
            old_profiles = self.get_iteration_points(check_points, self.__old_profiles[-1])

            for i in range(len(check_points)):

                for var in ["T", "P", "rho"]:

                    for j in [0, 1]:

                        abs_err = abs(old_profiles[j][i].get_variable(var) - curr_profiles[j][i].get_variable(var))
                        rel_err = abs_err / curr_profiles[j][i].get_variable(var)

                        if not (rel_err > ITERATION_TOLLS["rel"][var]):

                            return False

            return True

        return self.neglect_internal_heat_transfer

    def get_old_dh_int(self, is_annulus, depth):

        if self.neglect_internal_heat_transfer or (len(self.__old_profiles) == 0):

            return 0.

        try:

            p, rho = self.__old_profiles[-1][self.rw_geometry.i_annulus].get_iteration_value(depth)
            self.__tmp_point_annulus.set_variable("P", p)
            self.__tmp_point_annulus.set_variable("rho", rho)

            p, rho = self.__old_profiles[-1][self.rw_geometry.i_tubing].get_iteration_value(depth)
            self.__tmp_point_tubing.set_variable("P", p)
            self.__tmp_point_tubing.set_variable("rho", rho)

            q_int = self.rw_geometry.q_int(

                point_annulus=self.__tmp_point_annulus,
                point_tubing=self.__tmp_point_tubing

            )

            if is_annulus:

                return -q_int

            else:

                return q_int

        except:

            return 0.

    def get_profiles(self, position_list):

        profile_list = list()
        for profile in self.__old_profiles:
            profile_list.append(self.get_iteration_profile(position_list, profile))

        return profile_list