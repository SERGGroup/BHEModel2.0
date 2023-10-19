# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   IMPORT DATA                         -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "1 - No Geometry",
    "res"

)
filename = os.path.join(result_folder, "elapsed_times.npy")
times = np.load(filename)


# %%-------------------------------------   IMPORT DATA                         -------------------------------------> #
fig = plt.figure(figsize=(9, 4))
ax = fig.subplots(1, 1)
n_bins = 20

# Cumulative distributions.
n, bins, patches = ax.hist(

    times, n_bins, density=True, cumulative=False

)

plt.show()