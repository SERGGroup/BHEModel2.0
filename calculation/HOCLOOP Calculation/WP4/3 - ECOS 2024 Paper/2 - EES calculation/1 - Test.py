# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from EESConnect import EESConnector
from main_code import constants
import os.path

# %%------------   CALCULATIONS DATA                      -----------------------------------------------------------> #
t_sat = 10
t_max = 90
dt_he = 40
t_turb_in = 35
eta_comp = 0.8
eta_turb = 0.75
grad_T = 60
m_dot_well = 15


# %%------------   CALCULATE                              -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "2 - EES calculation",

)

EES_FILE = os.path.join(RES_FOLDER, "0 - EES Files", "base heat pump - python.EES")

with EESConnector(EES_FILE, ees_decimal_separator=",") as ees:

    try:
        result = ees.calculate([

            t_max, dt_he, t_turb_in,
            eta_comp, eta_turb, grad_T, m_dot_well

        ])
    except:

        result = [None, None, None]

    LCOH = result[0]
    q_net = result[1]
    w_net = result[2]
    t_max_reach = result[3]
    COP = q_net / w_net

print("Done!")
