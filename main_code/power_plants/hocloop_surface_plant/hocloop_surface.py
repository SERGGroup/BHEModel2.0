from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from EESConnect import EESConnector
import os.path


EES_FILE_DIR = os.path.join(os.path.dirname(__file__), "resources", "base heat pump - python.EES")

class HOCLOOPSurfaceHTHP:

    def __init__(

            self, hocloop_bhe: SimplifiedBHE,
            t_sat=10, t_max=50, dt_he=20,
            eta_comp=0.85, eta_turb=0.75,
            p_lim=20

    ):

        self.hocloop_bhe = hocloop_bhe

        self.t_sat = t_sat
        self.t_max = t_max
        self.dt_he = dt_he

        self.eta_comp = eta_comp
        self.eta_turb = eta_turb
        self.p_lim = p_lim

        self.w_net = None
        self.q_new = None
        self.x_out = None

    def update(self, time):

        self.hocloop_bhe.heating_section.time = time / 365
        self.hocloop_bhe.input_point.set_variable("T", self.t_sat)
        self.hocloop_bhe.input_point.set_variable("Q", 0.)
        self.hocloop_bhe.update()

        print("Well Updated!")

        with EESConnector(EES_FILE_DIR) as ees:

            try:
                result = ees.calculate([

                    self.t_sat,
                    self.t_max,
                    self.dt_he,
                    self.eta_comp,
                    self.eta_turb,
                    self.hocloop_bhe.output_point.get_variable("T"),
                    self.hocloop_bhe.output_point.get_variable("P") * 1000,
                    self.p_lim * 1000

                ])
            except:

                result = [None, None, None]

        self.w_net = result[0]
        self.q_new = result[1]
        self.x_out = result[2]


if __name__ == "__main__":

    from main_code.well_model.geometry_based_well_models.REELWEEL_model import REELWEELBHE, REELWELLGeometry, REELWELLRocksInfo
    from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
    from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
    import os

    #   Data from:
    #       "Deliverable 4.1 - Balmatt Case" (in "0 - Resources" Folder)
    #

    cas_id = 0.1617  # [m]

    depth = 2000  # [m]
    mass_flow = 8.8  # [kg/s]

    t_surf = 10  # [C]
    t_grad = 0.0325  # [C/m]
    t_res_top = 130  # [C]
    res_top = 3100  # [m]

    k_rock = 2.68  # [W/(m K)]
    c_rock = 0.93  # [kJ/(kg K)]
    rho_rock = 2600  # [kg/m^3]

    t_rock = t_res_top + t_grad * (depth - res_top)
    bhe_in = PlantThermoPoint(["Carbon Dioxide", "Ethylene"], [0.75, 0.25])
    bhe_in.set_variable("T", 20)
    bhe_in.set_variable("Q", 0)

    rocks_info = REELWELLRocksInfo(

        t_rocks=t_rock,
        k_rocks=k_rock,
        c_rocks=c_rock,
        rho_rocks=rho_rock,
        geo_gradient=t_grad,
        thermal_profile={

            0: 10,
            500: 25,
            900: 36,
            1600: 70,
            2200: 91,
            3100: 130,
            4500: 160

        }

    )

    hs_geometry = REELWELLGeometry(

        depth,
        tub_id=0.1,
        tub_od=0.13,
        cas_id=cas_id,
        cas_od=cas_id + 0.015,
        k_insulation=0.1,
        rocks_info=rocks_info,
        hot_in_tubing=True,
        max_back_time=3,
        alpha_old=0.5,
        neglect_internal_heat_transfer=True,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False

    )

    well = SimplifiedBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock,
        k_rocks=k_rock, c_rocks=c_rock, rho_rocks=rho_rock,
        t_surf=t_surf#, rw_geometry=hs_geometry, max_iteration=20

    )

    surf_hthp = HOCLOOPSurfaceHTHP(well, t_sat=-5)
    surf_hthp.update(7)

    print(well)