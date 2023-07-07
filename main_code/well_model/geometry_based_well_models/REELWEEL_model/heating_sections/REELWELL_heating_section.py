from main_code.well_model.simplified_well.heating_sections.abstract_class import AbstractHeatingSection
from main_code.support.other.simple_integrator import SimpleIntegrator
from scipy.integrate import RK45
from scipy.constants import g
import numpy as np


class REELWELLGeometry:

    def __init__(

            self, l_hor,
            tub_id=0.085, tub_od=0.16,
            cas_id=0.224, cas_od=0.244,
            k_insulation=0.15, ra_pipe=None,
            hot_in_tubing=False

    ):

        """

        REELWELL BHE well geometry.

        INPUT PARAMETERS:

            l_hor -> length of the horizontal section in [m]

            tub_id -> tubing internal diameter in [m]
            tub_od -> tubing external diameter in [m]

            cas_id -> casing internal diameter in [m]
            cas_od -> casing external diameter in [m]

            k_insulation -> thermal conductivity of the insulation between the tubing and the casing [W/(m*K)]
            ra_pipe -> roughness of the pipe surface [m]

        OTHER ATTRIBUTES:

            d_tub .......... diameter of the internal fluid tubing      [m]
            d_ann_int ...... internal diameter of the fluid annulus     [m]
            d_ann_out ...... external diameter of the fluid annulus     [m]
            d_h_ann ........ hydraulic diameter of the fluid annulus    [m]
            tkn_annulus .... thickness of the fluid annulus             [m]

            r_ins .......... thermal resistance of the insulation between tubing and annulus [(m*K)/W]
            parent_class ... reference to the REELWELLHeatingSection class
                             (needed for the ground heat exchange calculation)

        NOTE:

            - Default values for the geometry has been identified considering the REELWELL DualPipe
              dimensions

            - Natural Rubber has been considered in defining default value for thermal conductivity,
              this has to be clarified with REELWELL.

            - By default roughness is neglected (ra_pipe=None), correlation for smooth tubes are used.

        """

        self.tub_id = tub_id
        self.tub_od = tub_od

        self.cas_id = cas_id
        self.cas_od = cas_od

        self.l_hor = l_hor
        self.ra = ra_pipe
        self.k = k_insulation

        self.d_tub = self.tub_id
        self.d_ann_int = self.tub_od
        self.d_ann_out = self.cas_id

        self.d_h_ann = self.d_ann_out - self.d_ann_int
        self.tkn_annulus = (self.d_ann_out - self.d_ann_int) / 2

        self.r_ins = np.log(tub_od / tub_id) / (2 * np.pi * self.k)
        self.hot_in_tubing = hot_in_tubing
        self.parent_class = None

    def q_ground(self, point, depth=None):

        """

            This function returns the heat transfer between the annulus and the ground [kW/m]
            The implemented model is derived from:

                Y. Zhang, L. Pan, K. Pruess, S. Finsterle, Geothermics (2011), "A time-convolution approach for modeling heat
                exchange between a wellbore and surrounding formation"
                doi: 10.1016/j.geothermics.2011.08.003

        """

        if self.parent_class is not None:

            time = self.parent_class.time * 3.154e+7  # 3.154e+7 is a conversion factor: [year] -> [s]
            alpha = self.parent_class.main_BHE.alpha_rocks
            k = self.parent_class.main_BHE.k_rocks
            t_rock = self.get_t_rock(depth)
            dt = t_rock - point.get_variable("T")

            r0 = self.cas_od / 2
            td = alpha * time /  r0 ** 2

            if td < 2.8:

                f = k / r0 * ((np.pi * td) ** -0.5 + 0.5 - (td / np.pi) ** 0.5 / 4 + td / 8)

            else:

                gamma = 0.57722
                theta = (np.log(4 * td) - 2 * gamma)
                f = 2 * k / r0 * (1 / theta - gamma / theta ** 2)

            q_rel = f * dt / 1e3        # 1e3 is a conversion factor: [W/m^2] -> [kW/m^2]

            return q_rel * (2 * np.pi * r0)

        return 0.

    def q_int(self, point_annulus=None, point_tubing=None):

        """

            This function returns the heat transfer between the tubing and the annulus [kW/m]

            NOTE:

                !! The result is positive if heat is transferred from the annulus to the tubing
                !! point_annulus=None or point_tubing=None means that the heat transfer between the casing
                   and the tubing has been neglected

        """

        ###
        # OTHER NOTES:
        #
        #   - "1e3" is model conversion factor for "k" ([mW/(m*K)] - as returned by REFPROP -> [W/(m*K)])
        #   - "1e3" is model conversion factor for the returned value ([W/m] -> [kW/m])

        if point_tubing is None or point_annulus is None:

            return 0

        else:

            dt = point_annulus.get_variable("T") - point_tubing.get_variable("T")
            h_ann = self.nu(point_annulus, is_annulus=True) / self.d_h_ann * (point_annulus.get_variable("k") / 1e3)
            h_tub = self.nu(point_tubing, is_annulus=False) / self.d_tub * (point_tubing.get_variable("k") / 1e3)
            r = (1 / (self.d_tub * h_tub) + 1 / (self.d_ann_int * h_ann)) / np.pi + self.r_ins

            return dt / r / 1e3

    def nu(self, point, is_annulus):

        re = self.re(point, is_annulus)
        return 0.023 * np.power(re, 0.8) * np.power(point.get_variable("Pr"), 0.4)

    def re(self, point, is_annulus):

        # OTHER NOTES:
        #
        #   - "1e6" is a conversion factor for "mu" ([uPa/s] - as returned by REFPROP -> [Pa/s])

        if is_annulus:

            d_A = 4 / (np.pi * (self.d_ann_out + self.d_ann_int))

        else:

            d_A = 4 / (np.pi * self.d_tub)

        mu = point.get_variable("mu")

        if mu < 0:

            mu = 1e3

        return point.m_dot / mu * d_A * 1e6

    def f(self, point, is_annulus):

        # Friction Factor using Churchill Equation
        if self.ra is None:

            e_d = 0

        else:

            if is_annulus:

                e_d = self.ra / self.d_h_ann

            else:

                e_d = self.ra / self.d_tub

        re = self.re(point, is_annulus)
        A = (2.457 * np.log(1 / ((7. / re) ** 0.9 + 0.27 * e_d))) ** 16
        B = (37530. / re) ** 16
        f = 8. * ((8. / re) ** 12 + (1 / (A + B)) ** 1.5) ** (1. / 12.)

        return f

    def dp_dl(self, point, is_annulus):

        """

             This function returns the pressure derivative with length in the selected pipe [MPa/m]

        """

        ###
        # OTHER NOTES:
        #
        #   !! The fluid is considered to be locally incompressible !!
        #      (note the evaluation of the kinetic pressure)
        #
        #   - A is the flow area, dh is the hydraulic diameter
        #   - "1e6" is a conversion factor needed ([Pa] -> [MPa])
        #
        ###

        if is_annulus:

            A = np.pi * (self.d_ann_out - self.tkn_annulus) * self.tkn_annulus
            dh = self.d_h_ann

        else:

            A = np.pi / 4 * self.d_ann_out ** 2
            dh = self.d_tub

        rho = point.get_variable("rho")
        p_kinetic = point.m_dot ** 2 / (2 * rho * A ** 2)
        f = self.f(point, is_annulus)
        dp = - f / dh * p_kinetic  / 1e6

        return dp

    def dh_dl(self, point, is_annulus, other_point=None, depth=None):

        """

            This function returns the enthalpy derivative with length in the selected pipe [kJ/(kg * m)]

            NOTE:

                - If other_point is None means that the heat transfer between the casing
                  and the tubing has been neglected

        """

        if is_annulus:

            q_tot = self.q_ground(point, depth=depth) - self.q_int(point_annulus=point, point_tubing=other_point)

        else:

            q_tot = self.q_int(point_annulus=other_point, point_tubing=point)

        return q_tot / point.m_dot

    def get_t_rock(self, depth):

        t_rock = self.parent_class.main_BHE.t_rocks

        if depth is not None:

            grad = self.parent_class.main_BHE.geo_gradient
            t_rock -= (self.parent_class.main_BHE.dz_well - depth) * grad

        return t_rock

class REELWELLHeatingSection(AbstractHeatingSection):

    def __init__(

            self, main_BHE, reelwell_geometry: REELWELLGeometry,
            n_wells=1, neglect_internal_heat_transfer=True, time=1,
            integration_steps=None, integrate_temperature=False

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

        self.__tmp_ann = None
        self.__tmp_tub = None
        self.__tmp_point = None
        self.is_annulus = False

        self.integrator_list = list()
        self.integrators_profiler = list()
        self.integration_steps = integration_steps
        self.integrate_temperature = integrate_temperature

        self.neglect_internal_heat_transfer = neglect_internal_heat_transfer

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

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

        self.integrator_list = list()
        self.integrators_profiler = list()

        self.__tmp_ann = self.input_point.duplicate()
        self.__tmp_tub = self.input_point.duplicate()

        self.m_dot_well = self.main_BHE.m_dot / self.n_wells
        self.__tmp_ann.m_dot = self.m_dot_well
        self.__tmp_tub.m_dot = self.m_dot_well

        if self.neglect_internal_heat_transfer:

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

            self.integrator_list.append(RK45(

                self.integration_funct,
                pos_0, self.__convert_to_pt(x0), pos_end

            ))

        else:

            self.integrator_list.append(SimpleIntegrator(

                self.integration_funct,
                pos_0, self.__convert_to_pt(x0), pos_end,
                self.integration_steps

            ))

    def __integrate(self, pos_0, pos_end, x0, check_function=None):

        self.__append_integrator(pos_0, pos_end, x0)

        if len(self.integrator_list) > 0:

            __integrator = self.integrator_list[-1]
            self.integrators_profiler.append(list())

            while __integrator.status == 'running':

                __integrator.step()

                range_list = [__integrator.t_old, __integrator.t]
                range_list.sort()

                self.integrators_profiler[-1].append({

                    "range": range_list,
                    "dense_out": __integrator.dense_output(),
                    "error": (not __integrator.status == 'failed')

                })

                if check_function is not None:

                    if check_function(self.__convert_to_ph(self.integrator_list[-1].y)):

                        return None

            return self.__convert_to_ph(self.integrator_list[-1].y)

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
            dhdl = self.geom.dh_dl(self.__tmp_ann, self.is_annulus)

            if self.is_annulus ^ self.geom.hot_in_tubing:

                dxs = -dpdl, -dhdl

            else:

                dxs = dpdl, dhdl

        return dxs

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

    def get_heating_section_profile(self, position_list):

        t_list = np.full((2, len(position_list)), np.nan)
        p_list = np.full((2, len(position_list)), np.nan)

        if self.neglect_internal_heat_transfer:

            i_annulus = 0
            i_tubing = 1

            if not self.geom.hot_in_tubing:

                i_annulus = 1
                i_tubing = 0

            for i in range(len(position_list)):

                pos = position_list[i]

                p, h = self.get_heating_section_value(pos, i_annulus)
                self.__tmp_ann.set_variable("P", p)
                self.__tmp_ann.set_variable("H", h)

                p, h = self.get_heating_section_value(pos, i_tubing)
                self.__tmp_tub.set_variable("P", p)
                self.__tmp_tub.set_variable("H", h)

                t_list[0, i] = self.__tmp_ann.get_variable("T")
                t_list[1, i] = self.__tmp_tub.get_variable("T")
                p_list[0, i] = self.__tmp_ann.get_variable("P")
                p_list[1, i] = self.__tmp_tub.get_variable("P")

        return np.array(t_list), np.array(p_list)

    def get_heating_section_value(self, position, index):

        if len(self.integrators_profiler) > index:

            for step in self.integrators_profiler[index]:

                if step["range"][0] <= position <= step["range"][1]:

                    return self.__convert_to_ph(step["dense_out"](position))

            return self.__convert_to_ph(self.integrators_profiler[index][-1]["dense_out"](position))


class REELWELLInclinedHeatingSection(REELWELLHeatingSection):

    def __init__(

            self, main_BHE, reelwell_geometry: REELWELLGeometry,
            n_wells=1, neglect_internal_heat_transfer=True, time=1,
            integration_steps=None, integrate_temperature=False,
            inclination=90

    ):

        super().__init__(

            main_BHE, reelwell_geometry,
            n_wells=n_wells, neglect_internal_heat_transfer=neglect_internal_heat_transfer,
            time=time, integration_steps=integration_steps,
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
