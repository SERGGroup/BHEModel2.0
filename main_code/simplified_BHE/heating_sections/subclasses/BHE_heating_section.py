from main_code.simplified_BHE.heating_sections.abstract_class import AbstractHeatingSection
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
import numpy as np
import warnings


class BHEHeatingSection(AbstractHeatingSection):

    def __init__(

            self, main_BHE,
            DT_rocks=30, n_wells=10, d_well=0.5,
            tkn_annulus=0.05, save_profiles=False,
            consider_pressure_losses=True

    ):

        super().__init__(main_BHE)

        self.DT_rocks = DT_rocks
        self.n_wells = n_wells
        self.d_well = d_well
        self.tkn_annulus = tkn_annulus

        self.save_profiles = save_profiles
        self.consider_pressure_losses = consider_pressure_losses

        self.profiles = dict()

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def DT_rocks(self):

        return self.__DT_rocks

    @DT_rocks.setter
    def DT_rocks(self, DT_rocks_in):
        self.__DT_rocks = DT_rocks_in

    @property
    def n_wells(self):
        return self.__n_wells

    @n_wells.setter
    def n_wells(self, n_wells_in):
        self.__n_wells = n_wells_in

    @property
    def d_well(self):
        return self.__d_well

    @d_well.setter
    def d_well(self, d_well_in):
        self.__d_well = d_well_in

    @property
    def tkn_annulus(self):
        return self.__tkn_annulus

    @tkn_annulus.setter
    def tkn_annulus(self, tkn_annulus_in):
        self.__tkn_annulus = tkn_annulus_in

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_HS(self):

        """

            The Thermodynamic calculation proceeds as follows:

                - The HE section is divided in sections in which the fluid temperature increase by almost 1°C
                - The length of each section is calculated according to the heat transfer needs
                - the pressure losses are evaluated for the given length
                - the overall length and pressure loss is computed adding the lengths and the losses of each element

        """

        self.m_dot_well = self.main_BHE.m_dot / self.n_wells
        tmp_point_next = self.input_point.duplicate()
        tmp_point_old = self.input_point.duplicate()

        T_out = self.main_BHE.T_rocks - self.DT_rocks
        T_in = self.input_point.get_variable("T")

        self.l_HS = 0.
        self.dP_HS = 0.

        if T_out < T_in:

            # Impossible to extract work from such well, detailed analysis is unnecessary!

            self.output_point.set_variable("T", T_out)
            self.output_point.set_variable("P", self.input_point.get_variable("P"))

        else:

            n_sect = int(T_out - T_in)
            dT_sect = (T_out - T_in) / n_sect

            tmp_point_old.set_variable("P", self.input_point.get_variable("P"))
            tmp_point_old.set_variable("T", T_in)

            if self.save_profiles:

                self.__init_profiles(n_sect)
                self.__append_profile(0, 0, tmp_point_old)

            for i in range(n_sect):

                tmp_point_next.set_variable("P", tmp_point_old.get_variable("P"))
                tmp_point_next.set_variable("T", tmp_point_old.get_variable("T") + dT_sect)

                l_sect = self.__evaluate_section_length(tmp_point_old, tmp_point_next)
                dP_sect = self.__evaluate_section_pressure_losses(tmp_point_old, l_sect)

                self.l_HS += l_sect
                self.dP_HS += dP_sect

                tmp_point_old.set_variable("T", tmp_point_old.get_variable("T") + dT_sect)
                if self.consider_pressure_losses:
                    tmp_point_old.set_variable("P", self.input_point.get_variable("P") - self.dP_HS)
                else:
                    tmp_point_old.set_variable("P", self.input_point.get_variable("P"))

                if self.save_profiles:
                    self.__append_profile(i + 1, self.l_HS, tmp_point_old)

            self.output_point.set_variable("T", T_out)
            if self.consider_pressure_losses:
                self.output_point.set_variable("P", self.input_point.get_variable("P") - self.dP_HS)
            else:
                self.output_point.set_variable("P", self.input_point.get_variable("P"))

    def __evaluate_section_length(self, tmp_point_old, tmp_point_next):

        dT_rock_curr = self.main_BHE.T_rocks - tmp_point_old.get_variable("T")
        q_dot_need = self.m_dot_well * tmp_point_next.evaluate_variable_variation(tmp_point_old,"h")

        r_lin_ground = self.__calculate_ground_thermal_resistance()
        r_lin_fluid = self.__calculate_fluid_thermal_resistance(tmp_point_old)
        r_lin_tot = r_lin_ground + r_lin_fluid

        return (q_dot_need * 1E3 * r_lin_tot) / dT_rock_curr    # 1E3 is a conversion factor needed to convert kW in W

    def __calculate_ground_thermal_resistance(self):

        linear_resistance = 0.6161 / (self.main_BHE.geo_flux * (np.pi * self.d_well)) # TODO!! -> get right equation!!
        return linear_resistance

    def __calculate_fluid_thermal_resistance(self, tmp_point_old):

        NU_ext = self.__calculateNU()
        h_ext = NU_ext * (tmp_point_old.get_variable("k") / 1e3) / (2 * self.tkn_annulus)

        # Where:
        #
        #  - "1e3" is model conversion factor needed because REFPROP returns "k" in mW/(m*K) while we need it in W/(m*K)
        #  - "2 * self.tkn_annulus" is the hydraulic diameter
        #

        linear_resistance = 1. / (np.pi * self.d_well * h_ext)
        return linear_resistance

    def __calculateNU(self) -> float:

        # Nusselt correlation from Chengel "Heat Transfer model practical approach" pg.444. Values have been interpolated
        # from table 8-4 (Nusselt Number for fully developed -laminar- flow in an annulus with one surface isothermal
        # and the other -adiabatic-).
        #
        # Interpolation has been performed in an excel file that could be seen at this link:
        # https://firebasestorage.googleapis.com/v0/b/serg-group-repository.appspot.com/o/BHEModel%2Fannular%20flow%20correlation%20interpolation.xlsx?alt=media&token=5d644cf2-7aab-4878-9ead-506a38d49d19
        #
        # As can be seen, this correlation does not exactly match our needs because of the laminar flow an the adiabatic
        # wall requirements, further analysis have to be performed!!
        #
        # Tube roughness is taken into account following the approach described in Chengel "Heat Transfer model practical
        # approach" pg.443, hence considering the ratio between smooth and rough tube friction factor. Again this
        # method is not very accurate!

        x = (self.d_well - 2 * self.tkn_annulus) / self.d_well
        Nu_smooth = 4.0 + (0.85 - 0.008 * x) * x
        return Nu_smooth

    def __evaluate_section_pressure_losses(self, tmp_point_old, l_sect):

        f = self.__calculateFrictionFactor(tmp_point_old)
        rho = tmp_point_old.get_variable("rho")

        A = 2 * np.pi * (self.d_well - self.tkn_annulus) * self.tkn_annulus
        kynetic_pressure = self.m_dot_well ** 2 / (2 * rho * A ** 2)
        delta_P = l_sect / (2 * self.tkn_annulus) * f * kynetic_pressure / 1e6

        # Where:
        #
        #   - A is the flow area
        #   -
        #   - "1e6" is a conversion factor needed because Pressures are in in MPa while the resulting dP is in Pa
        #   - "2 * self.tkn_annulus" is the hydraulic diameter
        #

        return delta_P

    def __calculateFrictionFactor(self, tmp_point_old) -> float:

        # Friction Factor for smooth pipes calculated using Churchill Equation

        r_mean = (self.d_well - self.tkn_annulus) / 2
        Re = self.m_dot_well / (np.pi * r_mean * tmp_point_old.get_variable("mu"))

        A = (2.2113 * np.log(Re / 7.)) ** 16
        B = (37530. / Re) ** 16
        f = 8. * ((8. / Re) ** 12 + (1 / (A + B)) ** 1.5) ** (1. / 12.)

        return f

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def get_c_well(self):

        """

            Return the expected CAPEX for well drilling [€]

             Cost Correlation from:

                Adams et al, “Estimating the Geothermal Electricity Generation Potential of Sedimentary Basins Using
                genGEO (the generalizable GEOthermal techno-economic simulator)", ChemRxiv Prepr., 2021.

        """

        l_tot = self.l_HS * self.n_wells + self.main_BHE.dz_well
        c_well = (0.105 * l_tot ** 2 + 1776 * l_tot * self.d_well + 2.753 * pow(10, 5))

        if self.main_BHE.working_fluid == "CarbonDioxide":

            dc_CO2 = (133 + 256 * self.d_well) * l_tot

        else:

            dc_CO2 = 0.

        return c_well, dc_CO2

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def initial_guess(self):

        """
            Parameters to be optimized are:

                - DT_rocks
                - n_wells
                - d_well
                - tkn_annulus

        """
        return [50., 1, 0.2, 0.05]

    @property
    def optimization_bounds(self):

        """
            Parameters to be optimized are:

                - DT_rocks
                - n_wells
                - d_well
                - tkn_annulus

        """
        return Bounds([10., 1, 0.1, 0.005], [100., 50, 0.3, 0.05])

    def set_optimization_param(self, optimization_param):

        self.DT_rocks = optimization_param[0]
        self.n_wells = int(optimization_param[1])
        self.d_well = optimization_param[2]
        self.tkn_annulus = optimization_param[3]

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <---------------------   DISPLAY METHODS   ------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @staticmethod
    def __get_profile_parameter():

        return ["T", "P", "rho"]

    def __init_profiles(self, n_sect):

        self.profiles.update({"x": np.zeros((n_sect + 1, 1))})

        for key in self.__get_profile_parameter():

            self.profiles.update({key: np.zeros((n_sect + 1, 1))})

    def __append_profile(self, i, x, tmp_point):

        self.profiles["x"][i] = x

        for key in self.__get_profile_parameter():

            self.profiles[key][i] = tmp_point.get_variable(key)

    def plot_profiles(self):

        if not self.save_profiles:

            warnings.warn("!!IMPOSSIBLE TO DISPLAY THE PROFILE:\nThe 'save_profiles' flag was set to False when the calculation begins!!")
            return

        else:

            c_map = plt.cm.get_cmap("viridis")
            n_lines = len(self.__get_profile_parameter())
            fig, host = plt.subplots(figsize=(8, 5))  # (width, height) in inches

            par_list = [host]
            color_list = [c_map(0)]
            lns_list = list()

            for i in range(n_lines - 1):

                perc = (i + 1) / (n_lines - 1)
                par_list.append(host.twinx())
                color_list.append(c_map(perc))

            i = 0
            host.set_xlabel("Position [m]")
            for key in self.__get_profile_parameter():

                line, = par_list[i].plot(self.profiles["x"],self.profiles[key], color=color_list[i], label=key)
                par_list[i].set_ylabel("{} [{}]".format(key, self.input_point.get_unit(key)))
                par_list[i].yaxis.label.set_color(line.get_color())
                lns_list.append(line)
                i += 1

            if len(par_list) > 2:

                for i in range(2, len(par_list)):

                    par_list[i].spines['right'].set_position(('outward', 60 * (i - 1)))


            host.legend(handles=lns_list, loc='best')
            fig.tight_layout()
            plt.show()