# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os

if os.name == "nt":
    from main_code.constants import CALCULATION_FOLDER

else:
    CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations", "00 - Output")
water_folder = os.path.join(base_folder, "4 - Overall Optimization Result - Water", "2024-11-28 - Extended Optimization")
water_res = np.load(os.path.join(water_folder, 'optimized_array.npy'))
water_depth, water_grad = np.meshgrid(water_res[:, 0, 0, 0], water_res[0, :, 0, 1], indexing="ij")
water_points = np.vstack((water_depth.ravel(), water_grad.ravel())).T

co2_folder = os.path.join(base_folder, "3 - Overall Optimization Results", "2024-11-26 - Final Results")
co2_res = np.load(os.path.join(co2_folder, 'optimization_results.npy'))
co2_depth, co2_grad = np.meshgrid(co2_res[:, 0, 0, 0], co2_res[0, :, 0, 1], indexing="ij")
co2_points = np.vstack((co2_depth.ravel(), co2_grad.ravel())).T

depth_fine = np.linspace(np.nanmin(water_depth), np.nanmax(water_depth), 100)
grad_fine = np.linspace(np.nanmin(water_grad), np.nanmax(water_grad), 100)
x_fine, y_fine = np.meshgrid(depth_fine, grad_fine, indexing="ij")


# %%------------   PLOT OPTIMAL CONTOURS                  -----------------------------------------------------------> #
def frame_generation(frame):

    global fig, axs, water_res, co2_res, window_size, pbar

    for ax in axs:

        ax.clear()
        ax.set_xlabel("$T_{rocks}$ [°C]")

    omega = water_res[0, 0, frame, 2]

    if plot_compression:
        n = 0

    else:
        n = 1

    w_net = opt_res[:, :, frame, 4 + n * 3]
    m_dot = opt_res[:, :, frame, 5 + n * 3]

    t_sep_perc = opt_res[:, :, frame, 3]
    t_sep_perc_hthp = opt_res[:, :, frame, 6]

    if mean_in_omega:

        omega_i = frame + window_size
        omega_i += np.intc(np.linspace(-window_size, window_size, 2 * window_size + 1))

    else:

        omega_i = frame

    for i in range(len(x)):

        for j in range(len(y)):

            w_net[i, j] = np.nanmean(opt_res[i, j, omega_i, 4 + n * 3])
            m_dot[i, j] = np.nanmean(opt_res[i, j, omega_i, 5 + n * 3])

            t_sep_perc[i, j] = np.nanmean(opt_res[i, j, omega_i, 3])
            t_sep_perc_hthp[i, j] = np.nanmean(opt_res[i, j, omega_i, 6])

    if mean_in_space:

        w_net = gaussian_filter(griddata(points, w_net.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        m_dot = gaussian_filter(griddata(points, m_dot.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)

        t_sep_perc = gaussian_filter(griddata(points, t_sep_perc.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        t_sep_perc_hthp = gaussian_filter(griddata(points, t_sep_perc_hthp.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)

        x_plot = x_fine
        y_plot = y_fine

    else:

        x_plot = x
        y_plot = y

    if plot_optimal_depth:

        a = np.log(w_net + 1)
        b = np.log(1 / m_dot)
        x_min = (1 - omega) * a - omega * b

        optimal_depth = np.empty(len(y_fine))
        optimal_depth[:] = np.nan

        for j in range(len(y_fine)):
            i_min = np.where(x_min[:, j] == np.nanmax(x_min[:, j]))
            optimal_depth[j] = np.nanmean(x_fine[:, j][i_min])

    if plot_optim_params:
        axs[0].set_title("${T}_{sep\%}$ Compression [-]")
        axs[1].set_title("${T}_{sep\%}$ HTHP [-]")
        plot_elem = [t_sep_perc, t_sep_perc_hthp]

    else:
        axs[0].set_title("${W}_{net\ rel}$ [-]")
        axs[1].set_title("${m}_{ratio}$ [-]")
        plot_elem = [w_net, m_dot]

    for n in range(2):

        if levels[n] is None:
            contour = axs[n].contourf(y_plot, x_plot, plot_elem[n], extend='both', levels=21)
            levels[n] = contour.levels

        else:
            contour = axs[n].contourf(y_plot, x_plot, plot_elem[n], levels=levels[n], extend='both')

        if not bar_added[n]:

            fig.colorbar(contour, ax=axs[n])
            bar_added[n] = True

        if plot_optimal_depth:
            axs[n].plot(grad_fine, optimal_depth, lw=3, color='gold')

    fig.suptitle("Omega = {:.2f}".format(omega))

    try:
        pbar.update(1)
    except:
        pass
