# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.constants import CALCULATION_FOLDER
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
output_folder = os.path.join(base_folder, "00 - Output", "3 - Overall Optimization Results", "2024-11-26 - Final Results")

opt_res = np.load(os.path.join(output_folder, 'optimization_results.npy'))
x, y = np.meshgrid(opt_res[:, 0, 0, 0], opt_res[0, :, 0, 1], indexing="ij")
points = np.vstack((x.ravel(), y.ravel())).T

depth_fine = np.linspace(np.nanmin(x), np.nanmax(x), 100)
grad_fine = np.linspace(np.nanmin(y), np.nanmax(y), 100)
x_fine, y_fine = np.meshgrid(depth_fine, grad_fine, indexing="ij")


# %%------------   PLOT OPTIMAL CONTOURS                  -----------------------------------------------------------> #
plot_optim_params = True
generate_gif = True

mean_in_omega = False
window_size = 0

mean_in_space = True
interp_method = "cubic"
sigma_smoot = 15

plot_optimal_depth = False
add_w_net_zero_contour = True
add_ihx_zero_contour = True
ihx_zero_contour_sp = 0.25

def frame_generation(frame):

    global fig, axs, opt_res, window_size, pbar

    for ax in axs:

        ax.clear()

        ax.set_xlabel("Geothermal Gradient [°C/km]")
        ax.set_ylabel("Depth [km]")
        ax.invert_yaxis()

    omega = opt_res[0, 0, frame, 2]

    T_max_perc = opt_res[:, :, 0, 3]
    sep_perc = opt_res[:, :, 0, 4]
    w_net = opt_res[:, :, 0, 5]
    m_dot = opt_res[:, :, 0, 6]
    ihx_power = opt_res[:, :, 0, 9]

    if mean_in_omega:

        omega_i = frame + window_size
        omega_i += np.intc(np.linspace(-window_size, window_size, 2 * window_size + 1))

    else:

        omega_i = frame


    for i in range(len(x)):

        for j in range(len(y)):

            T_max_perc[i, j] = np.nanmean(opt_res[i, j, omega_i, 3])
            sep_perc[i, j] = np.nanmean(opt_res[i, j, omega_i, 4])
            w_net[i, j] = np.nanmean(opt_res[i, j, omega_i, 5])
            m_dot[i, j] = np.nanmean(opt_res[i, j, omega_i, 6])
            ihx_power[i, j] = np.nanmean(opt_res[i, j, omega_i, 9])

    if mean_in_space:

        T_max_perc = gaussian_filter(griddata(points, T_max_perc.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        sep_perc = gaussian_filter(griddata(points, sep_perc.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        w_net = gaussian_filter(griddata(points, w_net.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        m_dot = gaussian_filter(griddata(points, m_dot.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)
        ihx_power = gaussian_filter(griddata(points, ihx_power.ravel(), (x_fine, y_fine), method=interp_method), sigma=sigma_smoot)

        x_plot = x_fine
        y_plot = y_fine

    else:

        x_plot = x
        y_plot = y

    ihx_power = (ihx_power - np.nanmin(ihx_power)) / (np.nanmax(ihx_power) - np.nanmin(ihx_power))
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
        axs[0].set_title("${T}_{SG\%}$ [-]")
        axs[1].set_title("${sep}_{\%}$ [-]")
        plot_elem = [T_max_perc, sep_perc]

    else:
        axs[0].set_title("${W}_{net\ rel}$ [-]")
        axs[1].set_title("${m}_{ratio}$ [-]")
        plot_elem = [w_net, m_dot]

    for n in range(2):

        if levels[n] is None:
            contour = axs[n].contourf(y_plot, x_plot, plot_elem[n], extend='both')
            levels[n] = contour.levels

        else:
            contour = axs[n].contourf(y_plot, x_plot, plot_elem[n], levels=levels[n], extend='both')

        if not bar_added[n]:

            fig.colorbar(contour, ax=axs[n])
            bar_added[n] = True

        if add_w_net_zero_contour:
            axs[n].contour(y_plot, x_plot, w_net, levels=[0], linewidths=3, colors='white')

        if add_ihx_zero_contour:
            axs[n].contour(y_plot, x_plot, ihx_power, levels=[ihx_zero_contour_sp], linewidths=3, colors='tab:orange')

        if plot_optimal_depth:
            axs[n].plot(grad_fine, optimal_depth, lw=3, color='gold')

    fig.suptitle("Omega = {:.2f}".format(omega))

    try:
        pbar.update(1)
    except:
        pass

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
bar_added = [False, False]

if not generate_gif:

    levels = [21, 21]
    frame_generation(68)
    plt.show()

else:

    fps = 20
    remove_n_omega = 50
    levels = [np.linspace(0, 1, 21), np.linspace(0, 1, 21)]

    if mean_in_omega:
        n_frames = opt_res.shape[2] - (2 * window_size + 1) - remove_n_omega

    else:
        n_frames = opt_res.shape[2] - 1 - remove_n_omega
        window_size = 0

    pbar = tqdm(total=n_frames+1)
    ani = FuncAnimation(fig, frame_generation, frames=n_frames, repeat=True)

    if plot_optim_params:
        ani.save(os.path.join(output_folder, f'optimization results - params.gif'), writer='pillow', fps=fps)

    else:
        ani.save(os.path.join(output_folder, f'optimization results.gif'), writer='pillow', fps=fps)

    plt.close(fig)
    pbar.close()


# %%------------   PLOT X_min OVER DEPTH                  -----------------------------------------------------------> #
grad_plot = np.array([20, 40, 60, 80])
plot_animation = False

mean_in_space = True
interp_method = "cubic"
sigma_smoot = 10

def opt_frame_generation(frame):

    global fig, ax, opt_res, window_size, pbar, x_plot, y_plot, legend_added, grad_plot

    omega = opt_res[0, 0, frame, 2]

    w_net = opt_res[:, :, 0, 5]
    m_dot = opt_res[:, :, 0, 6]

    for i in range(len(x)):

        for j in range(len(y)):

            w_net[i, j] = np.nanmean(opt_res[i, j, frame, 5])
            m_dot[i, j] = np.nanmean(opt_res[i, j, frame, 6])

    w_net = gaussian_filter(griddata(points, w_net.ravel(), (x_plot, y_plot), method=interp_method), sigma=sigma_smoot)
    m_dot = gaussian_filter(griddata(points, m_dot.ravel(), (x_plot, y_plot), method=interp_method), sigma=sigma_smoot)

    a = np.log(w_net + 1)
    b = np.log(1 / m_dot)
    x_min = (1 - omega) * a - omega * b

    ax.clear()

    for k in range(len(grad_plot)):

        j_curr = np.where(y_plot == grad_plot[k])
        ax.plot(x_plot[j_curr]/1e3, x_min[j_curr], label="{}".format(grad_plot[k]))

    if not legend_added:
       ax.legend()
       legend_added = True

    ax.set_ylabel("x_min [-]")
    ax.set_xlabel("Depth [km]")
    fig.suptitle("Omega = {:.2f}".format(omega))

    try:
        pbar.update(1)
    except:
        pass

fig, ax = plt.subplots(1, 1, figsize=(12, 5))
legend_added = False

grad_plot_final = np.concatenate((grad_fine, grad_plot))
grad_plot_final = np.sort(np.unique(grad_plot_final))
x_plot, y_plot = np.meshgrid(depth_fine, grad_plot_final, indexing="ij")

if plot_animation:

    fps = 20
    remove_n_omega = 50
    n_frames = opt_res.shape[2] - 1 - remove_n_omega

    pbar = tqdm(total=n_frames+1)

    ani = FuncAnimation(fig, opt_frame_generation, frames=n_frames, repeat=True)
    ani.save(os.path.join(output_folder, f'x_min over grad results.gif'), writer='pillow', fps=fps)

    plt.close(fig)
    pbar.close()

else:

    opt_frame_generation(68)
    plt.show()