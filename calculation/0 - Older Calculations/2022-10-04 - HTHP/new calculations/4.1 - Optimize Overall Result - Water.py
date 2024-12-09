# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os

if os.name == "nt":
    from main_code.constants import CALCULATION_FOLDER

else:
    CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
output_folder = os.path.join(base_folder, "00 - Output", "4 - Overall Optimization Result - Water", "2024-11-28 - Extended Optimization")

n_omega = 149
omega_list = np.linspace(0, 1, n_omega)

output_array = np.load(os.path.join(output_folder, 'output_array.npy'))
depth_list = np.unique(output_array[:, :, :, 0])
grad_list = np.unique(output_array[:, :, :, 1])
t_sep_list = np.unique(output_array[:, :, :, 2])

n_sep_fine = 100

depth_mesh, grad_mesh = np.meshgrid(depth_list, grad_list, indexing="ij")
new_shape = (depth_mesh.shape[0], depth_mesh.shape[1], n_omega, 9)
optimized_array = np.empty(new_shape)
optimized_array[:] = np.nan


# %%------------   FIND OPTIMUM                           -----------------------------------------------------------> #
pbar = tqdm(total=len(depth_list) * len(grad_list))
cmap = get_cmap("viridis")

ignore_plotting = False
plot_only_result = True
plot_every = 40
interp_kind = 'linear'

a = 0
sigma_smooth = 5
plot_original = True
original_result = optimized_array[0, 0, :, :]

for i, depth in enumerate(depth_list):

    optimized_array[i, :, :, 0] = depth

    for j, grad in enumerate(grad_list):

        optimized_array[i, j, :, 1] = grad
        optimized_array[i, j, :, 2] = omega_list

        plot_result = ((i * len(grad_list) + j) % plot_every == 0) and (not ignore_plotting)
        if plot_result and (not plot_only_result):

            fig, axs = plt.subplots(1, 2)
            axs[0].set_title("COMPRESSION SYSTEM")
            axs[1].set_title("HTHP")

            plt.suptitle(f"Depth: {depth}m, Grad: {grad}°C/km")

        for n in range(2):

            w_net_curr = output_array[i, j, :, 4 + 2 * n]
            m_ratio_curr = output_array[i, j, :, 3 + 2 * n]

            acptb_curr = np.where(np.logical_and(

                m_ratio_curr > 0,
                w_net_curr > -1

            ))

            t_sep_curr = t_sep_list[acptb_curr]
            w_net_curr = w_net_curr[acptb_curr]
            m_ratio_curr = m_ratio_curr[acptb_curr]

            t_sep_fine = np.linspace(np.nanmin(t_sep_curr), np.nanmax(t_sep_curr), n_sep_fine)
            interp_w_net = interp1d(t_sep_curr, w_net_curr, kind=interp_kind, fill_value="extrapolate")
            interp_m_ratio = interp1d(t_sep_curr, m_ratio_curr, kind=interp_kind, fill_value="extrapolate")

            w_net_curr = interp_w_net(t_sep_fine)
            m_ratio_curr = interp_m_ratio(t_sep_fine)

            a_curr = np.log(w_net_curr + 1)
            b_curr = np.log(1 / m_ratio_curr)

            for k, omega in enumerate(omega_list):

                x_min_curr = (1 - omega) * a_curr - omega * b_curr
                i_min_curr = np.where(x_min_curr == np.nanmax(x_min_curr))
                optimized_array[i, j, k, 3 + 3 * n] = np.nanmean(t_sep_fine[i_min_curr])
                optimized_array[i, j, k, 4 + 3 * n] = np.nanmean(w_net_curr[i_min_curr])
                optimized_array[i, j, k, 5 + 3 * n] = np.nanmean(m_ratio_curr[i_min_curr])

                if plot_result and (not plot_only_result):

                    axs[n].plot(t_sep_fine, x_min_curr, color=cmap(omega))
                    axs[n].plot(

                        t_sep_curr[i_min_curr], x_min_curr[i_min_curr],
                        marker='*', markersize=15, color='#FFD700'

                    )

            original_result[:, 3 + 3 * n: 6 + 3 * n] = optimized_array[i, j, :, 3 + 3 * n: 6 + 3 * n]

            optimized_array[i, j, :, 3 + 3 * n] = gaussian_filter1d(optimized_array[i, j, :, 3 + 3 * n], sigma=sigma_smooth)
            optimized_array[i, j, :, 4 + 3 * n] = interp_w_net(optimized_array[i, j, :, 3 + 3 * n])
            optimized_array[i, j, :, 5 + 3 * n] = interp_m_ratio(optimized_array[i, j, :, 3 + 3 * n])

        if plot_result:

            if not plot_only_result:
                plt.show()

            fig, ax = plt.subplots(1, 1)

            modifier = 1

            if a == 1:
                modifier = -1

            line_comp, = ax.plot(omega_list, modifier*optimized_array[i, j, :, 3 + a], label="COMPRESSION SYSTEM")
            line_hp, = ax.plot(omega_list, modifier*optimized_array[i, j, :, 6 + a], label="HTHP")

            if plot_original:

                ax.plot(omega_list, modifier*original_result[:, 3 + a], "--", color=line_comp.get_color())
                ax.plot(omega_list, modifier*original_result[:, 6 + a], "--", color=line_hp.get_color())

            ax.set_ylabel("$T_{sep\ curr}$ [-]")

            if a == 1:
                ax.set_ylabel("$W_{net}$ [-]")

            elif a == 2:
                ax.set_ylabel("$m_{dot} [-]")
                ax.set_yscale("log")

            ax.set_xlabel("$\Omega$ [-]")
            ax.legend()

            plt.suptitle(f"Depth: {depth}m, Grad: {grad}°C/km")
            plt.show()

        np.save(os.path.join(output_folder, "optimized_array.npy"), optimized_array)
        pbar.update(1)

pbar.close()
