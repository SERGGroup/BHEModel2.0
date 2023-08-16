from main_code.well_model.simplified_well.heating_sections.abstract_class import AbstractHeatingSection
from main_code.well_model.geometry_based_well_models.REELWEEL_model.geometry import REELWELLGeometry
from main_code.support.other.integration_profiler import IntegrationProfiler
from main_code.support.other.simple_integrator import SimpleIntegrator
from scipy.integrate import RK45
from scipy.constants import g
from copy import deepcopy
import numpy as np


class REELWELLHeatingSection(AbstractHeatingSection):

    def __init__(

            self, main_BHE, reelwell_geometry: REELWELLGeometry,
            n_wells=1, time=1, integration_steps=None,
            integrate_temperature=False

    ):

        super().__init__(main_BHE)

        self.n_wells = n_wells
        self.geom = reelwell_geometry
        self.geom.parent_class = self

        self.m_dot_well = 0.
        self.dt_rocks = 0.

        self.n_bisect = 20
        self.bisect_toll = 1e-5
        self.time = time

        self.is_annulus = False

        self.__tmp_ann = None
        self.__tmp_tub = None
        self.__tmp_point = None
        self.__old_profiles = None

        self.integrators_profiler = list()
        self.integration_steps = integration_steps
        self.integrate_temperature = integrate_temperature

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def reset_old_profilers(self):

        self.__old_profiles = list()

    def update_HS(self):

        """

            The Thermodynamic calculation proceeds as follows:

                - If internal heat transfer is neglected, a Runge-Kutta integration is performed first along the tubing
                  then in the casing. The RK function returns the derivative of pressure and enthalpy through the pipe.

                - If internal heat transfer is considered, the two branches of the pipe are solved simultaneously. The
                  outlet conditions are iterated in order to equalize the conditions at the end of the tubing and at
                  the beginning of the casing

                - If hot_in_tubing is true, the fluid will pass first into the annulus and then into the tubing.
                  The solution order will be switched by changing the __is_annulus property:

                        self.__is_annulus = self.hot_in_tubing -> True if self.hot_in_tubing == True
                        self.__is_annulus = not self.hot_in_tubing -> False if self.hot_in_tubing == True

        """

        if len(self.integrators_profiler) > 0:
            self.__old_profiles.append(deepcopy(self.integrators_profiler))

        self.integrators_profiler = list()

        self.__tmp_ann = self.input_point.duplicate()
        self.__tmp_tub = self.input_point.duplicate()

        self.m_dot_well = self.main_BHE.m_dot / self.n_wells
        self.__tmp_ann.m_dot = self.m_dot_well
        self.__tmp_tub.m_dot = self.m_dot_well

        if self.solve_sequentially:

            # if the internal heat transfer is neglected or the main_BHE class is the REELWEELBHE class
            # (recognized by the fact that it has an attribute called rw_geometry), the integration
            # proceeds in a standard way

            self.is_annulus = self.geom.hot_in_tubing
            x0 = (self.input_point.get_variable("P"), self.input_point.get_variable("h"))
            x_mid = self.__integrate(pos_0=0, pos_end=self.geom.l_hor, x0=x0)

            self.is_annulus = not self.geom.hot_in_tubing
            res = self.__integrate(pos_0=self.geom.l_hor, pos_end=0, x0=x_mid)

        else:

            res = self.__bisect_output()

        self.output_point.set_variable("P", res[0])
        self.output_point.set_variable("h", res[1])

    def __convert_to_pt(self, x0):

        if self.integrate_temperature:

            if self.__tmp_point is None:
                self.__tmp_point = self.input_point.duplicate()

            self.__tmp_point.set_variable("P", x0[0])
            self.__tmp_point.set_variable("h", x0[1])
            x0_1 = self.__tmp_point.get_variable("T")

            if len(x0) > 2:

                self.__tmp_point.set_variable("P", x0[2])
                self.__tmp_point.set_variable("h", x0[3])
                x0_3 = self.__tmp_point.get_variable("T")

                x0 = (x0[0], x0_1, x0[2], x0_3)

            else:

                x0 = (x0[0], x0_1)

        return x0

    def __convert_to_ph(self, x0):

        if self.integrate_temperature:

            if self.__tmp_point is None:
                self.__tmp_point = self.input_point.duplicate()

            self.__tmp_point.set_variable("P", x0[0])
            self.__tmp_point.set_variable("T", x0[1])
            x0_1 = self.__tmp_point.get_variable("h")

            if len(x0) > 2:

                self.__tmp_point.set_variable("P", x0[2])
                self.__tmp_point.set_variable("T", x0[3])
                x0_3 = self.__tmp_point.get_variable("h")

                x0 = (x0[0], x0_1, x0[2], x0_3)

            else:

                x0 = (x0[0], x0_1)

        return x0

    def __convert_derivatives(self, x0, dxs):

        if self.integrate_temperature:

            if self.__tmp_point is None:
                self.__tmp_point = self.input_point.duplicate()

            self.__tmp_point.set_variable("P", x0[0])
            self.__tmp_point.set_variable("H", x0[1])

            dhdp = self.__tmp_point.get_derivative("H", "P", "T")
            dhdt = self.__tmp_point.get_derivative("H", "T", "P")

            dx_1 = (dxs[1] - dhdp * dxs[0]) / dhdt

            if len(x0) > 2:

                self.__tmp_point.set_variable("P", x0[2])
                self.__tmp_point.set_variable("H", x0[3])

                dhdp = self.__tmp_point.get_derivative("H", "P", "T")
                dhdt = self.__tmp_point.get_derivative("H", "P", "T")

                dx_3 = (dxs[3] - dhdp * dxs[2]) / dhdt

                dxs = (dxs[0], dx_1, dxs[2], dx_3)

            else:

                dxs = (dxs[0], dx_1)

        return dxs

    def __append_integrator(self, pos_0, pos_end, x0):

        if self.integration_steps is None:

            self.integrators_profiler.append(
                IntegrationProfiler(
                    RK45(
                        self.integration_funct,
                        pos_0, self.__convert_to_pt(x0), pos_end
                    )
                )
            )

        else:

            self.integrators_profiler.append(
                IntegrationProfiler(
                    SimpleIntegrator(
                        self.integration_funct,
                        pos_0, self.__convert_to_pt(x0), pos_end,
                        self.integration_steps
                    )
                )
            )

    def __integrate(self, pos_0, pos_end, x0, check_function=None):

        self.__append_integrator(pos_0, pos_end, x0)

        if len(self.integrators_profiler) > 0:

            __integrator = self.integrators_profiler[-1]

            while __integrator.status == 'running':

                __integrator.step()

                if check_function is not None:

                    if check_function(self.__convert_to_ph(self.integrators_profiler[-1].y)):

                        return None

            return self.__convert_to_ph(self.integrators_profiler[-1].y)

        return None

    def __bisect_output(self):

        t_interval, p_interval = self.get_bisect_intervals()
        res = []

        for i in range(self.n_bisect):

            t_guess = np.mean(t_interval)
            p_guess = np.mean(p_interval)

            self.__tmp_ann.set_variable("T", t_guess)
            self.__tmp_ann.set_variable("P", p_guess)

            x0 = (

                self.__tmp_ann.get_variable("P"),
                self.__tmp_ann.get_variable("h"),
                self.input_point.get_variable("P"),
                self.input_point.get_variable("h")

            )

            self.integrators_profiler = list()
            res = self.__integrate(pos_0=0, pos_end=self.geom.l_hor, x0=x0, check_function=self.bisect_check)

            if res is None:

                t_interval[0] = t_guess
                p_interval[0] = p_guess

            else:

                p_ann, p_tub, t_ann, t_tub = self.__get_check_values(res)
                p_err = abs(p_ann - p_tub) / p_tub
                t_err = abs(t_ann - t_tub) / t_tub

                if  p_err < self.bisect_toll and t_err < self.bisect_toll:

                    break

                else:


                    i = (1, 0)[t_ann < t_tub]
                    t_interval[i] = t_guess

                    i = (1, 0)[p_ann < p_tub]
                    p_interval[i] = p_guess

        return res

    def __get_check_values(self, res):

        p_ann = res[0]
        p_tub = res[2]

        self.__tmp_ann.set_variable("P", p_ann)
        self.__tmp_ann.set_variable("h", res[1])
        self.__tmp_tub.set_variable("P", p_tub)
        self.__tmp_tub.set_variable("h", res[3])

        t_ann = self.__tmp_ann.get_variable("T")
        t_tub = self.__tmp_tub.get_variable("T")

        return p_ann, p_tub, t_ann, t_tub

    def get_bisect_intervals(self):

        t_interval = [self.input_point.get_variable("T"), self.main_BHE.t_rocks]
        p_interval = [0., self.input_point.get_variable("P")]

        return t_interval, p_interval

    def bisect_check(self, y):

        p_ann, p_tub, t_ann, t_tub = self.__get_check_values(y)

        if t_ann < t_tub and p_ann < p_tub:

            return True

        return False

    def integration_funct(self, t, x):

        """

            Runge-Kutta Integration Function:

                The function returns the derivative of pressure and enthalpy in the pipe.
                Depending on the number of elements of x different calculation are performed:

                    - If x is a list with 4 elements both fluid condition in both the annulus and the tubing are solved
                      simultaneously. The first two elements are considered to be pressure and enthalpy in the annulus,

                    - If x has only 2 elements fluid condition are solved one after the other in the annulus and in the
                      tubing. The inputs are pressure and enthalpy of the tubing or of the annulus depending on the
                      "is_annulus" arguments

        """

        # NOTE:
        #
        #   The function returns negative values for the derivatives in the annulus because the integration proceeds
        #   upstream in that section. If hot_in_tubing is true all the changes are switched because the fluid has
        #   changed direction in the pipes.
        #
        #   Note that this has been implemented with the XOR operator (^).
        #
        #       "self.__is_annulus ^ self.hot_in_tubing" is true if one between __is_annulus or hot_in_tubing is True
        #       (but not if they both are)
        #

        return self.__convert_derivatives(x, self.base_integration_funct(t, self.__convert_to_ph(x)))

    def base_integration_funct(self, t, x):

        """

            Runge-Kutta Integration Function:

                The function returns the derivative of pressure and enthalpy in the pipe.
                Depending on the number of elements of x different calculation are performed:

                    - If x is a list with 4 elements both fluid condition in both the annulus and the tubing are solved
                      simultaneously. The first two elements are considered to be pressure and enthalpy in the annulus,

                    - If x has only 2 elements fluid condition are solved one after the other in the annulus and in the
                      tubing. The inputs are pressure and enthalpy of the tubing or of the annulus depending on the
                      "is_annulus" arguments

        """

        # NOTE:
        #
        #   The function returns negative values for the derivatives in the annulus because the integration proceeds
        #   upstream in that section. If hot_in_tubing is true all the changes are switched because the fluid has
        #   changed direction in the pipes.
        #
        #   Note that this has been implemented with the XOR operator (^).
        #
        #       "self.__is_annulus ^ self.hot_in_tubing" is true if one between __is_annulus or hot_in_tubing is True
        #       (but not if they both are)
        #

        if len(x) == 4:

            self.__tmp_ann.set_variable("P", x[0])
            self.__tmp_ann.set_variable("h", x[1])
            self.__tmp_tub.set_variable("P", x[2])
            self.__tmp_tub.set_variable("h", x[3])

            dpdl_ann = self.geom.dp_dl(self.__tmp_ann, is_annulus=True)
            dpdl_tub = self.geom.dp_dl(self.__tmp_tub, is_annulus=False)

            dhdl_ann = self.geom.dh_dl(self.__tmp_ann, is_annulus=True, other_point=self.__tmp_tub)
            dhdl_tub = self.geom.dh_dl(self.__tmp_tub, is_annulus=False, other_point=self.__tmp_ann)

            if self.geom.hot_in_tubing:

                dxs = dpdl_ann, dhdl_ann, -dpdl_tub, -dhdl_tub

            else:

                dxs = -dpdl_ann, -dhdl_ann, dpdl_tub, dhdl_tub

        else:

            self.__tmp_ann.set_variable("P", x[0])
            self.__tmp_ann.set_variable("h", x[1])

            dpdl = self.geom.dp_dl(self.__tmp_ann, self.is_annulus)
            dh_ext = self.geom.dh_dl(self.__tmp_ann, self.is_annulus)
            dh_int = self.get_old_dh_int(is_annulus=self.is_annulus, depth=t)

            if self.is_annulus ^ self.geom.hot_in_tubing:

                dhdl = dh_ext + dh_int
                dxs = -dpdl, -dhdl

            else:

                dhdl = dh_ext + dh_int
                dxs = dpdl, dhdl

        return dxs

    @property
    def solve_sequentially(self):

        return self.geom.neglect_internal_heat_transfer or (not self.main_BHE.neglect_internal_heat_transfer)

    def get_old_dh_int(self, is_annulus, depth):

        if self.geom.neglect_internal_heat_transfer or (len(self.__old_profiles) == 0):

            return 0.

        try:

            if self.integrate_temperature:

                q_int = self.geom.get_old_q_int(

                    is_annulus=is_annulus, depth=depth,
                    old_profiles=self.__old_profiles,
                    first_var="P", second_var="T"

                )

            else:

                q_int = self.geom.get_old_q_int(

                    is_annulus=is_annulus, depth=depth,
                    old_profiles=self.__old_profiles,
                    first_var="P", second_var="T"

                )

            if is_annulus:

                return -q_int / self.__tmp_ann.m_dot

            else:

                return q_int / self.__tmp_tub.m_dot

        except:

            return 0.

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def get_c_well(self):
        # TODO
        return 0.

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def initial_guess(self):
        return []

    @property
    def optimization_bounds(self):
        return []

    def set_optimization_param(self, optimization_param):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <---------------------   DISPLAY METHODS   ------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def get_heating_section_profile(self, position_list, profile=None):

        if profile is None:
            profile = self.integrators_profiler

        t_list = np.full((2, len(position_list)), np.nan)
        p_list = np.full((2, len(position_list)), np.nan)
        rho_list = np.full((2, len(position_list)), np.nan)
        h_list = np.full((2, len(position_list)), np.nan)

        if self.solve_sequentially:

            i_annulus = self.geom.i_annulus
            i_tubing = self.geom.i_tubing

            for i in range(len(position_list)):

                pos = position_list[i]

                if self.integrate_temperature:

                    p, t = profile[i_annulus].get_iteration_value(pos)
                    self.__tmp_ann.set_variable("P", p)
                    self.__tmp_ann.set_variable("T", t)

                    p, t = profile[i_tubing].get_iteration_value(pos)
                    self.__tmp_tub.set_variable("P", p)
                    self.__tmp_tub.set_variable("T", t)

                else:

                    p, h = profile[i_annulus].get_iteration_value(pos)
                    self.__tmp_ann.set_variable("P", p)
                    self.__tmp_ann.set_variable("H", h)

                    p, h = profile[i_tubing].get_iteration_value(pos)
                    self.__tmp_tub.set_variable("P", p)
                    self.__tmp_tub.set_variable("H", h)

                t_list[0, i] = self.__tmp_ann.get_variable("T")
                t_list[1, i] = self.__tmp_tub.get_variable("T")
                p_list[0, i] = self.__tmp_ann.get_variable("P")
                p_list[1, i] = self.__tmp_tub.get_variable("P")

                rho_list[0, i] = self.__tmp_ann.get_variable("rho")
                rho_list[1, i] = self.__tmp_tub.get_variable("rho")
                h_list[0, i] = self.__tmp_ann.get_variable("H")
                h_list[1, i] = self.__tmp_tub.get_variable("H")

        return np.array(t_list), np.array(p_list), np.array(rho_list), np.array(h_list)

    def additional_setup_data(self, data_frame: dict):

        internal_he = "ignored"
        integration_type = "RK45"

        if self.integration_steps is not None:

            integration_type = "simple - {} steps".format(self.integration_steps)

        if not self.geom.neglect_internal_heat_transfer:

            internal_he = "considered"

        data_frame["Calculation Options"].update({

            "heating section": {"value": "REELWELLHeatingSection", "unit": None},
            "hs - integration": {"value": integration_type, "unit": None},
            "hs - internal heat exc": {"value": internal_he, "unit": None},

        })

        data_frame = self.geom.additional_setup_data(data_frame, is_well=False)

        return data_frame

    def get_profiles(self, position_list, get_index=None):

        if get_index is None:

            profile_list = list()
            for profile in self.__old_profiles:
                profile_list.append(self.get_heating_section_profile(position_list, profile))

            return profile_list

        else:

            profile = self.__old_profiles[get_index]
            return self.get_heating_section_profile(position_list, profile)


class REELWELLInclinedHeatingSection(REELWELLHeatingSection):

    def __init__(

            self, main_BHE, reelwell_geometry: REELWELLGeometry,
            n_wells=1, time=1, integration_steps=None,
            integrate_temperature=False,
            inclination=90

    ):

        super().__init__(

            main_BHE, reelwell_geometry,
            n_wells=n_wells, time=time, integration_steps=integration_steps,
            integrate_temperature=integrate_temperature


        )

        self.inclination = inclination
        self.__tmp_ann = None
        self.__tmp_tub = None

    @property
    def inclination(self):

        return self.__alpha / np.pi * 180

    @inclination.setter
    def inclination(self, inclination):

        self.__alpha = inclination / 180 * np.pi

    def base_integration_funct(self, t, x):

        """

            Runge-Kutta Integration Function:

                The function returns the derivative of pressure and enthalpy in the pipe.
                Depending on the number of elements of x different calculation are performed:

                    - If x is a list with 4 elements both fluid condition in both the annulus and the tubing are solved
                      simultaneously. The first two elements are considered to be pressure and enthalpy in the annulus,

                    - If x has only 2 elements fluid condition are solved one after the other in the annulus and in the
                      tubing. The inputs are pressure and enthalpy of the tubing or of the annulus depending on the
                      "is_annulus" arguments

                    - Pressure and enthalpy increase due to gravitational potential change is considered as the current
                      integrator represent a vertical well

        """

        # NOTE:
        #
        #   The function returns negative values for the derivatives in the annulus because the integration proceeds
        #   upstream in that section. If hot_in_tubing is true all the changes are switched because the fluid has
        #   changed direction in the pipes. The gravitational contribution is positive in any case as it increase with
        #   depth
        #
        #   Note that this has been implemented with the XOR operator (^).
        #
        #       "self.__is_annulus ^ self.hot_in_tubing" is true if one between __is_annulus or hot_in_tubing is True
        #       (but not if they both are)
        #

        depth = t * np.sin(self.__alpha)
        dh_grav = g * np.sin(self.__alpha) / 1e3

        if self.__tmp_ann is None:

            self.__tmp_ann = self.input_point.duplicate()
            self.__tmp_tub = self.input_point.duplicate()

        if len(x) == 4:

            self.__tmp_ann.set_variable("P", x[0])
            self.__tmp_ann.set_variable("h", x[1])
            self.__tmp_tub.set_variable("P", x[2])
            self.__tmp_tub.set_variable("h", x[3])

            dp_grav_ann = dh_grav * self.__tmp_ann.get_variable("rho") / 1e3
            dp_grav_tub = dh_grav * self.__tmp_tub.get_variable("rho") / 1e3

            dpdl_ann = self.geom.dp_dl(self.__tmp_ann, is_annulus=True)
            dpdl_tub = self.geom.dp_dl(self.__tmp_tub, is_annulus=False)

            dhdl_ann = self.geom.dh_dl(self.__tmp_ann, is_annulus=True, other_point=self.__tmp_tub, depth=depth)
            dhdl_tub = self.geom.dh_dl(self.__tmp_tub, is_annulus=False, other_point=self.__tmp_ann, depth=depth)

            if self.geom.hot_in_tubing:

                return dpdl_ann + dp_grav_ann, dhdl_ann + dh_grav, -dpdl_tub + dp_grav_tub, -dhdl_tub + dh_grav

            return -dpdl_ann + dp_grav_ann, -dhdl_ann + dh_grav, dpdl_tub + dp_grav_tub, dhdl_tub + dh_grav

        else:

            self.__tmp_ann.set_variable("P", x[0])
            self.__tmp_ann.set_variable("h", x[1])
            dp_grav = g * np.sin(self.__alpha) * self.__tmp_ann.get_variable("rho") / 1e6

            dpdl = self.geom.dp_dl(self.__tmp_ann, self.is_annulus)
            dhdl = self.geom.dh_dl(self.__tmp_ann, self.is_annulus, depth=depth)

            # if self.__tmp_ann.get_variable("T") < 0:
            #
            #     print("T")

            if self.is_annulus ^ self.geom.hot_in_tubing:

                return -dpdl + dp_grav, -dhdl + dh_grav

            return dpdl + dp_grav, dhdl + dh_grav

    def get_bisect_intervals(self):

        # TODO

        t_interval = [self.input_point.get_variable("T"), self.main_BHE.t_rocks]
        p_interval = [0., self.input_point.get_variable("P")]

        return t_interval, p_interval

    def additional_setup_data(self, data_frame: dict):

        data_frame = super().additional_setup_data(data_frame)

        data_frame["Calculation Options"].update({

            "heating section": {"value": "REELWELLHeatingSection", "unit": None}

        })

        data_frame["well geometry"].update({

            "inclination": {"value": self.inclination, "unit": "°"}

        })

        return data_frame
