from main_code.well_model.ground_models import evaluate_ground_f
import numpy as np


ITERATION_TOLLS = {

    "rel": {

        "T": 0.001,
        "P": 0.001,
        "rho": 0.001

    }

}

class REELWELLRocksInfo:

    def __init__(

            self, t_rocks = None,
            k_rocks = None, alpha_rocks = None,
            c_rocks = None, rho_rocks = None,
            geo_gradient = None, thermal_profile = None

    ):

        self.__t_rocks = t_rocks
        self.__k_rocks = k_rocks
        self.__c_rocks = c_rocks
        self.__rho_rocks = rho_rocks

        if (alpha_rocks is None) and (not None in [k_rocks, c_rocks, rho_rocks]):

            self.__alpha_rocks = k_rocks / (rho_rocks * c_rocks * 1e3)

        else:

            self.__alpha_rocks = alpha_rocks

        self.__geo_gradient = geo_gradient
        self.__thermal_profile = thermal_profile

        self.main_class = None

    @property
    def t_rocks(self):

        if self.__t_rocks is not None:
            return self.__t_rocks

        if not self.main_class.parent_class is None:
            return self.main_class.parent_class.main_BHE.t_rocks

        return None

    @property
    def k_rocks(self):

        if self.__k_rocks is not None:
            return self.__k_rocks

        if not self.main_class.parent_class is None:
            return self.main_class.parent_class.main_BHE.k_rocks

        return None

    @property
    def alpha_rocks(self):

        if self.__alpha_rocks is not None:
            return self.__alpha_rocks

        if not self.main_class.parent_class is None:
            return self.main_class.parent_class.main_BHE.alpha_rocks

        return None

    @property
    def geo_gradient(self):

        if self.__geo_gradient is not None:
            return self.__geo_gradient

        if not self.main_class.parent_class is None:
            return self.main_class.parent_class.main_BHE.geo_gradient

        return None

    @t_rocks.setter
    def t_rocks(self, t_rocks_in):

        self.__t_rocks = t_rocks_in

    @k_rocks.setter
    def k_rocks(self, k_rocks_in):

        self.__k_rocks = k_rocks_in

    @alpha_rocks.setter
    def alpha_rocks(self, alpha_rocks_in):

        self.__alpha_rocks = alpha_rocks_in

    @geo_gradient.setter
    def geo_gradient(self, geo_gradient_in):

        self.__geo_gradient = geo_gradient_in

    def get_t_rock(self, depth = None):

        t_rock = self.t_rocks

        if depth is not None:

            if self.__thermal_profile is not None:

                depths = list(self.__thermal_profile.keys())
                temperatures = list(self.__thermal_profile.values())

                if depth > depths[-1]:

                    grad = self.geo_gradient
                    t_rock = grad * (depth - depths[-1]) + temperatures[-1]

                else:

                    for i in range(len(depths)):

                        if depths[i] > depth:

                            grad = (temperatures[i] - temperatures[i - 1]) / (depths[i] - depths[i - 1])
                            t_rock = grad * (depth - depths[i - 1]) + temperatures[i - 1]
                            break

            else:

                grad = self.geo_gradient
                t_rock -= (self.main_class.parent_class.main_BHE.dz_well - depth) * grad

        return t_rock


class REELWELLGeometry:

    def __init__(

            self, l_hor,
            tub_id=0.085, tub_od=0.16,
            cas_id=0.224, cas_od=0.244,
            k_insulation=0.15, ra_pipe=None,
            max_back_time=5, alpha_old = 0.5,
            rocks_info=None, hot_in_tubing=False, neglect_internal_heat_transfer=True,
            ignore_tubing_pressure_losses=False,
            ignore_annulus_pressure_losses=False

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

            CALCULATION PARAMETERS:

                1) Internal heat transfer calculation options
                max_back_time -> maximum number of previous profiles conserved in memory
                alpha_old -> exponential impact of old profiles on the calculation

                2) Nusselt modification options
                nu_modifier -> modifier which increase or decrease the nusselt
                               number to simulate differences with heat transfer predictions


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

        self.d_hyd_ann = self.d_ann_out - self.d_ann_int
        self.tkn_annulus = (self.d_ann_out - self.d_ann_int) / 2

        self.a_flow_tub = np.pi / 4 * self.d_tub ** 2
        self.a_flow_ann = np.pi / 4 * (self.d_ann_out ** 2 - self.d_ann_int ** 2)

        self.r_ins = np.log(tub_od / tub_id) / (2 * np.pi * self.k)

        self.hot_in_tubing = hot_in_tubing
        self.max_back_time = max_back_time
        self.alpha_old = alpha_old

        if rocks_info is None:
            self.rocks_info = REELWELLRocksInfo()
        else:
            self.rocks_info = rocks_info
        self.rocks_info.main_class = self

        self.neglect_internal_heat_transfer = neglect_internal_heat_transfer
        self.ignore_tubing_pressure_losses = ignore_tubing_pressure_losses
        self.ignore_annulus_pressure_losses = ignore_annulus_pressure_losses

        self.parent_class = None
        self.__tmp_point_ann = None
        self.__tmp_point_tub = None

    def q_ground(self, point, depth=None):

        if self.parent_class is not None:

            k = self.rocks_info.k_rocks
            time = self.parent_class.time * 3.154e+7  # 3.154e+7 is a conversion factor: [year] -> [s]
            alpha = self.rocks_info.alpha_rocks

            r0 = self.cas_od / 2
            td = alpha * time /  r0 ** 2

            f = evaluate_ground_f(td)

            r_rocks = r0 / (k * f)

            k_fluid = point.get_variable("k") / 1e3
            r_fluid = self.d_hyd_ann / (self.nu(point, is_annulus=True) * k_fluid)

            r_tot = r_rocks + r_fluid

            t_rock = self.rocks_info.get_t_rock(depth)
            dt = t_rock - point.get_variable("T")
            q_rel = dt / r_tot / 1e3  # 1e3 is a conversion factor: [W/m^2] -> [kW/m^2]

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

            h_ann = self.nu(point_annulus, is_annulus=True) * (point_annulus.get_variable("k") / 1e3) / self.d_hyd_ann
            h_tub = self.nu(point_tubing, is_annulus=False) * (point_tubing.get_variable("k") / 1e3) / self.d_tub

            r = (1 / (self.d_tub * h_tub) + 1 / (self.d_ann_int * h_ann)) / np.pi + self.r_ins

            return dt / r / 1e3

    def nu(self, point, is_annulus):

        re = self.re(point, is_annulus)
        pr = point.get_variable("Pr")
        re_pow = np.sign(re) * np.power(np.abs(re), 0.8)

        if type(pr) == float:

            pr_pow = np.sign(pr) * np.power(np.abs(pr), 0.4)
            return 0.023 * re_pow * pr_pow

        else:
            return 0.023 * re_pow * np.power(0.7, 0.4)

    def re(self, point, is_annulus):

        # OTHER NOTES:
        #
        #   - "1e6" is a conversion factor for "mu" ([uPa/s] - as returned by REFPROP -> [Pa/s])

        if is_annulus:

            d_A = self.d_hyd_ann / self.a_flow_ann

        else:

            d_A = self.d_tub / self.a_flow_tub

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

                e_d = self.ra / self.d_hyd_ann

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

            if self.ignore_annulus_pressure_losses:

                return 0.

            A = self.a_flow_ann
            d_hyd = self.d_hyd_ann

        else:

            if self.ignore_tubing_pressure_losses:

                return 0.

            A = self.a_flow_tub
            d_hyd = self.d_tub

        rho = point.get_variable("rho")
        p_kinetic = point.m_dot ** 2 / (2 * rho * A ** 2)
        f = self.f(point, is_annulus)
        dp = - f / d_hyd * p_kinetic  / 1e6

        return dp

    def dh_dl(self, point,  is_annulus, other_point=None, depth=None):

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

    def get_old_q_int(self, is_annulus, depth, old_profiles, first_var="P", second_var="rho"):

        if self.neglect_internal_heat_transfer or (len(old_profiles) == 0):

            return 0.

        try:

            if self.__tmp_point_ann is None:
                self.__tmp_point_ann = self.parent_class.main_BHE.input_point.duplicate()

            if self.__tmp_point_tub is None:
                self.__tmp_point_tub = self.parent_class.main_BHE.input_point.duplicate()

            q_int = 0.
            sum_mean = 0.
            for i in range(self.max_back_time):

                if i < len(old_profiles):

                    var_1, var_2 = old_profiles[-(i+1)][self.i_annulus].get_iteration_value(depth)
                    self.__tmp_point_ann.set_variable(first_var, var_1)
                    self.__tmp_point_ann.set_variable(second_var, var_2)

                    var_1, var_2 = old_profiles[-(i+1)][self.i_tubing].get_iteration_value(depth)
                    self.__tmp_point_tub.set_variable(first_var, var_1)
                    self.__tmp_point_tub.set_variable(second_var, var_2)

                    q_int_curr = self.q_int(

                        point_annulus=self.__tmp_point_ann,
                        point_tubing=self.__tmp_point_tub

                    )

                else:

                    q_int_curr = 0.

                sum_mean += self.alpha_old ** i
                q_int += q_int_curr * self.alpha_old ** i

            q_int = q_int / sum_mean

            if is_annulus:

                return -q_int

            else:

                return q_int

        except:

            return 0.

    def additional_setup_data(self, data_frame: dict, is_well=True):

        if is_well:

            key = "well geometry"

        else:

            key = "Heating Section Data"

        data_frame[key].update({

            "tub_id": {"value": self.tub_id, "unit": "m"},
            "tub_od": {"value": self.tub_od, "unit": "m"},
            "cas_id": {"value": self.cas_id, "unit": "m"},
            "cas_od": {"value": self.cas_od, "unit": "m"},
            "k_insulation": {"value": self.k, "unit": "W/(m*K)"},
            "d_tub": {"value": self.d_tub, "unit": "m"},
            "d_ann_int": {"value": self.d_ann_int, "unit": "m"},
            "d_ann_out": {"value": self.d_ann_out, "unit": "m"},
            "tkn_annulus": {"value": self.tkn_annulus, "unit": "m"},
            "d_h_ann": {"value": self.d_hyd_ann, "unit": "m"},
            "r_ins": {"value": self.r_ins, "unit": "(m*K)/W"}

        })

        return data_frame

    def fluid_speed(self, is_annulus, point):

        if is_annulus:

            A = self.a_flow_ann

        else:

            A = self.a_flow_tub

        return point.m_dot / (point.get_variable("rho") * A)

    @property
    def i_annulus(self):

        if not self.hot_in_tubing:
            return 1
        return 0

    @property
    def i_tubing(self):

        if not self.hot_in_tubing:
            return 0
        return 1
