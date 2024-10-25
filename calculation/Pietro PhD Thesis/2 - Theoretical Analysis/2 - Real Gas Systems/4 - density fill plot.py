# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER
from REFPROPConnector import ThermodynamicPoint
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import colors
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   CALCULATIONS OPTIONS                -------------------------------------> #
n_points = 500

t_rels = np.linspace(0.5, 2, n_points)
p_rels = np.logspace(-2, 2, n_points)

t_tp_mesh, p_tp_mesh = np.meshgrid(t_rels, p_rels, indexing='ij')


# %%-------------------------------------   INIT FLUID                          -------------------------------------> #
fluid = "Carbon Dioxide"

tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
tp.set_variable("T", tp.RPHandler.TC)
tp.set_variable("P", tp.RPHandler.PC)
cp0_crit = tp.get_variable("CP0")
rho_crit = tp.get_variable("rho")

t_crit = tp.RPHandler.TC
p_crit = tp.RPHandler.PC


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
res_tp_mesh = np.empty(p_tp_mesh.shape)
p_sat_arr = np.empty((len(t_rels), 1))
res_tp_mesh[:] = np.nan
p_sat_arr[:] = np.nan

pbar = tqdm(desc="Calculating Points", total=n_points ** 2)

for i in range(len(t_rels)):

    t_curr = t_crit * t_rels[i]

    if t_rels[i] < 1:

        tp.set_variable("T", t_curr)
        tp.set_variable("Q", 0)
        p_sat_arr[i, 0] = tp.get_variable("P") / p_crit

    for k in range(len(p_rels)):

        p_curr = p_rels[k] * p_crit
        tp.set_variable("T", t_curr)
        tp.set_variable("P", p_curr)

        res_tp_mesh[i, k] = tp.get_variable("rho") / rho_crit

        pbar.update(1)

pbar.close()


# %%-------------------------------------   PLOT CONTOUR                        -------------------------------------> #
fig, ax = plt.subplots(1, 1, figsize=(7, 4), dpi=600)

n_level = 100
n_ticks_level = 5
contour_alpha = 0.6

ax.set_yscale("log")
ax.set_xlabel("$T_{red}$ [-]")
ax.set_ylabel("$p_{red}$ [-]")
ax.set_xlim(t_rels.min(), t_rels.max())
ax.set_ylim(p_rels.min(), p_rels.max())

# Customizing the secondary x-axis
ax2 = ax.twiny()
ax2.set_xlabel("$T_{CO_2}$ [°C]")
ax2.set_xlim(t_rels.min()*t_crit - 273.15, t_rels.max()*t_crit-273.15)

ax3 = ax.twinx()
ax3.set_ylabel("$P_{CO_2}$ [bar]")
ax3.set_ylim(p_rels.min()*p_crit/1e5, p_rels.max()*p_crit/1e5)
ax3.set_yscale("log")

levels = np.logspace(-1, 0.5, n_level) * rho_crit
ticks_levels =  [50, 100, 225, 500, 1000]

# Create a colormap with alpha
cmap = plt.cm.viridis  # Choose your colormap
cmap_alpha = cmap(np.linspace(1, 0, cmap.N))
inv_colormap = colors.ListedColormap(cmap_alpha)

ax.plot(t_rels, p_sat_arr[:, 0], color="black", linewidth=2, zorder=2)
cs = ax.contourf(

    t_tp_mesh, p_tp_mesh, res_tp_mesh * rho_crit,
    levels=levels, zorder=1, extend='both',
    cmap = inv_colormap

)
# Add a semi-transparent white rectangle to simulate alpha effect
rect = Rectangle((0, 0), 1, 1, transform=ax.transAxes, color='white', alpha=1-contour_alpha)
ax.add_patch(rect)

# Manually adjust the position of the main axes
plt.subplots_adjust(left=0.1, right=0.75, top=0.85, bottom=0.15)

# Add colorbar as a separate axis
cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(cs, cax=cbar_ax, alpha=contour_alpha)
cbar.set_ticks(ticks_levels)
cbar.ax.set_ylabel("$\\rho_{CO_2}$ [kg/$m^3$]")

plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filename = "rho_behaviour_contour.png"
filepath = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "3.general output", filename

)
fig.savefig(filepath)
