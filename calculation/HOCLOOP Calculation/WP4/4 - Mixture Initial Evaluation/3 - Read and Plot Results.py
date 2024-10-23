# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from tkinter.filedialog import askopenfilename
from pandas import DataFrame, ExcelWriter
from main_code import constants
from datetime import datetime
import numpy as np
import os


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "4 - Mixture Initial Evaluation",
    "res"

)

filepath = askopenfilename(initialdir=RES_FOLDER)
results = np.load(filepath)


# %%------------   GENERATE EXCEL FILE                    -----------------------------------------------------------> #
now = datetime.now()
now_str = "{}-{:02}-{:02} - {:02}.{:02}".format(now.year, now.month, now.day, now.hour, now.minute)
file_path = os.path.join(RES_FOLDER, "excel-results", "performance-results-{}.xlsx".format(now_str))
df = DataFrame(results)
df.columns = [

    "Grad_Geo [°C/km]", "Flow Rate [kg/s]",
    "concentration [%mass]", "time [days]",
    "beta well [-]", "dT well [°C]",
    "dh well [kJ/kg]", "dex well [kJ/kg]",
    "dex bottom [kJ/kg]", "well power [kW]",
    "well exergy increase [kW]",
    "LCOH [€/kJ]", "LCOex [€/kJ]"

]

df.to_excel(file_path, index=False, startrow=1, startcol=1)
