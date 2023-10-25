# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   IMPORT DATA                         -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "0.support", "res"

)
filename = os.path.join(result_folder, "elapsed_times.npy")
times = np.load(filename, allow_pickle=True)


# %%-------------------------------------   IMPORT DATA                         -------------------------------------> #
fig = plt.figure(figsize=(9, 4))
ax = fig.subplots(1, 1)
n_bins = 500

# Cumulative distributions.
n, bins, patches = ax.hist(

    times, n_bins, density=True, cumulative=False

)
# plt.xlim((300, 350))
plt.show()