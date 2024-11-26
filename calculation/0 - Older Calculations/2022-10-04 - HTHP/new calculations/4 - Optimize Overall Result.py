# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.ndimage import gaussian_filter, uniform_filter
from main_code.constants import CALCULATION_FOLDER
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_laplace
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import warnings
import os


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
# CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
output_folder = os.path.join(base_folder, "00 - Output", "3 - Overall Optimization Results")
gif_folder = os.path.join(base_folder, "00 - Output", "4 - Control GIFs")

depth_list = np.load(os.path.join(output_folder, 'depth_list.npy'))
grad_list = np.load(os.path.join(output_folder, 'grad_list.npy'))
t_sg_perc_list = np.load(os.path.join(output_folder, 't_sg_perc_list.npy'))
sep_perc_list = np.load(os.path.join(output_folder, 'sep_perc_list.npy'))

x, y = np.meshgrid(t_sg_perc_list, sep_perc_list, indexing="ij")
xy_points = np.vstack((x.ravel(), y.ravel())).T

w_rel_max = np.load(os.path.join(output_folder, 'w_rel_max.npy'))
m_dot_max = np.load(os.path.join(output_folder, 'm_dot_max.npy'))
p_max_list = np.load(os.path.join(output_folder, 'p_max_list.npy'))
t_ihx_list = np.load(os.path.join(output_folder, 't_ihx_list.npy'))
ihx_power_list = np.load(os.path.join(output_folder, 'ihx_power_list.npy'))


# %%------------   FIND OPTIMUM                           -----------------------------------------------------------> #
sigma_val = 3
window_size = 5
window_size_filter = 5  # Should be an odd number
poly_order = 2   # Polynomial order

gaussian_filter_before = True
gaussian_filter_after = False
use_savgol = True

fps = 50
s_len = 3

n_samples = 20
dev_perc = 0.001
n_frames = fps * s_len

optimal_results = np.empty((

    len(depth_list), len(grad_list),
    n_frames, 10

))
optimal_results[:] = np.nan
plot_gif = False
plot_pareto = True

pbar = tqdm(total=len(depth_list) * len(grad_list))

for i, depth in enumerate(depth_list):

    for j, grad in enumerate(grad_list):

        with warnings.catch_warnings():

            warnings.simplefilter("ignore", category=RuntimeWarning)
            a = np.log(w_rel_max[i, j, :, :] + 1)
            b = np.log(1 / m_dot_max[i, j, :, :])

        if gaussian_filter_before:

            a = gaussian_filter(a, sigma=sigma_val)
            b = gaussian_filter(b, sigma=sigma_val)

        if plot_gif:
            fig, axs = plt.subplots(1, 3, figsize=(12, 5))

        if plot_pareto:

            fig, ax = plt.subplots(1, 1, figsize=(6, 5))

            m_ratio_flat = m_dot_max[i, j, :, :].flatten()
            w_rel_flat = w_rel_max[i, j, :, :].flatten()
            t_sg_perc_flat = x[:, :].flatten()
            sep_perc_flat = y[:, :].flatten()

            otm_points = np.vstack((m_ratio_flat, w_rel_flat)).T
            pareto_points = {

                "X": list(),
                "Y": list(),

                "m_ratio": list(),
                "w_net": list()

            }

            median_calculation = []

            for m, curr_point in enumerate(otm_points):

                is_pareto = True
                for k, other_point in enumerate(otm_points):

                    if all(other_point >= curr_point) and any(other_point > curr_point):
                        is_pareto = False
                        break

                if is_pareto:

                    pareto_points["X"].append(t_sg_perc_flat[m])
                    pareto_points["Y"].append(sep_perc_flat[m])
                    pareto_points["m_ratio"].append(m_ratio_flat[m])
                    pareto_points["w_net"].append(w_rel_flat[m])

                    if 0.1 < sep_perc_flat[m] < 0.8 and t_sg_perc_flat[m] > t_sg_perc_flat[0]:

                        median_calculation.append(t_sg_perc_flat[m])

            x_mean = np.mean(median_calculation)
            ax.scatter(pareto_points["X"], pareto_points["Y"], alpha=0.5)

            x_array = np.concatenate((

                np.linspace(np.min(x), x_mean, 10),
                np.ones(30)*x_mean,
                np.linspace(x_mean, np.max(x), 10)

            ))

            y_array = np.concatenate((

                np.ones(10) * np.max(y),
                np.linspace(np.min(x), np.max(x), 30),
                np.ones(10) * np.min(y),

            ))

            ax.plot(x_array, y_array, linewidth=3, color="tab:green", alpha=0.5)

            w_rel_opt = griddata(xy_points, w_rel_max[i, j].flatten(), (x_array, y_array), method='cubic')
            m_dot_opt = griddata(xy_points, m_dot_max[i, j].flatten(), (x_array, y_array), method='cubic')
            p_max_opt = griddata(xy_points, p_max_list[i, j].flatten(), (x_array, y_array), method='cubic')
            t_ihx_opt = griddata(xy_points, t_ihx_list[i, j].flatten(), (x_array, y_array), method='cubic')
            ihx_power_opt = griddata(xy_points, ihx_power_list[i, j].flatten(), (x_array, y_array), method='cubic')

            a_opt = np.log(w_rel_opt + 1)
            b_opt = np.log(1 / m_dot_opt)

            for k, omega in enumerate(np.linspace(0, 1, n_frames - 1)):

                optimal_results[i, j, k, 0] = depth
                optimal_results[i, j, k, 1] = grad
                optimal_results[i, j, k, 2] = omega

                x_min = (1 - omega) * a_opt - omega * b_opt
                x_min_ovr = np.nanmin(x_min)

                k_min = np.where(x_min == x_min_ovr)
                x_star = np.nanmean(x_array[k_min])
                y_star = np.nanmean(y_array[k_min])

                optimal_results[i, j, k, 3] = x_star
                optimal_results[i, j, k, 4] = y_star

                optimal_results[i, j, k, 5] = np.nanmean(w_rel_opt[k_min])
                optimal_results[i, j, k, 6] = np.nanmean(m_dot_opt[k_min])
                optimal_results[i, j, k, 7] = np.nanmean(p_max_opt[k_min])
                optimal_results[i, j, k, 8] = np.nanmean(t_ihx_opt[k_min])
                optimal_results[i, j, k, 9] = np.nanmean(ihx_power_opt[k_min])

        def frame_generation(frame):

            global a, b, depth, grad
            global n_frames, optimal_results

            omega = (frame + 1) / n_frames

            x_min = (1 - omega) * a - omega * b
            dx = savgol_filter(x_min, window_size_filter, poly_order, deriv=1, axis=1, delta=1)
            dy = savgol_filter(x_min, window_size_filter, poly_order, deriv=1, axis=0, delta=1)
            dtot = dx * dx + dy * dy

            if gaussian_filter_after:
                dtot = gaussian_filter(dtot, sigma=sigma_val)

            min_der = np.nanmin(dtot)
            optimal_results[i, j, frame, 0] = depth
            optimal_results[i, j, frame, 1] = grad
            optimal_results[i, j, frame, 2] = omega

            if min_der < 1e-5:

                i_min = np.where(dtot == min_der)

                x_star = np.nanmean(x[i_min])
                y_star = np.nanmean(y[i_min])

                optimal_results[i, j, frame, 3] = x_star
                optimal_results[i, j, frame, 4] = y_star

                optimal_results[i, j, frame, 5] = w_rel_max[i, j][i_min]
                optimal_results[i, j, frame, 6] = m_dot_max[i, j][i_min]
                optimal_results[i, j, frame, 7] = p_max_list[i, j][i_min]
                optimal_results[i, j, frame, 8] = t_ihx_list[i, j][i_min]
                optimal_results[i, j, frame, 9] = ihx_power_list[i, j][i_min]

            if plot_gif:

                axs[0].clear()
                axs[0].contourf(x, y, x_min, levels=15, cmap='viridis')
                axs[0].set_title("x_min")

                axs[1].clear()
                axs[1].contourf(x, y, np.log10(dtot), levels=15, cmap='viridis')
                axs[1].set_title("dtot")

                log_image = gaussian_laplace(x_min, sigma=1)
                axs[2].clear()
                axs[2].contourf(x, y, log_image, levels=15, cmap='viridis')
                axs[2].contour(x, y, log_image, levels=[0], colors='white', linewidths=2)
                axs[2].set_title("Gaussian Laplace")

                fig.suptitle("{}m - {}C - Omega {:.2f}".format(depth, grad, omega))

                if min_der < 1e-5:

                    for ax in axs:

                        ax.plot(x_star, y_star, marker='*', markersize=15, color='#FFD700')

        if plot_gif:

            ani = FuncAnimation(fig, frame_generation, frames=n_frames-1, repeat=True)
            ani.save(os.path.join(gif_folder, f'{depth}m - {grad}C-m.gif'), writer='pillow', fps=fps)
            np.save(os.path.join(output_folder, 'optimization_results.npy'), optimal_results)

            plt.close(fig)

        if plot_pareto:

            # for k in range(n_frames - 1):
            #
            #     frame_generation(k)

            ax.plot(optimal_results[i, j, :, 3], optimal_results[i, j, :, 4], alpha=0.5, linewidth=2, color="tab:orange")
            plt.savefig(os.path.join(gif_folder, f'{depth}m - {grad}C-m.png'), format='png', dpi=300)
            plt.close(fig)

            np.save(os.path.join(output_folder, 'optimization_results.npy'), optimal_results)

        pbar.update(1)

pbar.close()
