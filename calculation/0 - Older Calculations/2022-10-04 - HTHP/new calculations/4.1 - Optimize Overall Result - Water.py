# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
# from main_code.constants import CALCULATION_FOLDER
from scipy.interpolate import interp1d
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os

CALCULATION_FOLDER = "/Users/PietroUngar/PycharmProjects/BHEModel2.0/calculation"


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
base_folder = os.path.join(CALCULATION_FOLDER, "0 - Older Calculations", "2022-10-04 - HTHP", "new calculations")
output_folder = os.path.join(base_folder, "00 - Output", "4 - Overall Optimization Result - Water", "2024-11-28 - Extended Optimization")

n_omega = 50
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

        for k, omega in enumerate(omega_list):

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

                interp_w_net = interp1d(t_sep_curr, w_net_curr, kind='linear')
                interp_m_ratio = interp1d(t_sep_curr, m_ratio_curr, kind='linear')

                t_sep_fine = np.linspace(np.nanmin(t_sep_curr), np.nanmax(t_sep_curr), n_sep_fine)
                w_net_curr = interp_w_net(t_sep_fine)
                m_ratio_curr = interp_m_ratio(t_sep_fine)

                a_curr = np.log(w_net_curr + 1)
                b_curr = np.log(1 / m_ratio_curr)
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

        if plot_result:

            if not plot_only_result:
                plt.show()

            fig, ax = plt.subplots(1, 1)

            a = 2
            ax.plot(omega_list, optimized_array[i, j, :, 3 + a], label="COMPRESSION SYSTEM")
            ax.plot(omega_list, optimized_array[i, j, :, 6 + a], label="HTHP")

            ax.set_xlabel("$\Omega$ [-]")
            ax.set_ylabel("$T_{sep\ curr}$ [-]")
            ax.legend()

            plt.suptitle(f"Depth: {depth}m, Grad: {grad}°C/km")
            plt.show()

        np.save(os.path.join(output_folder, "optimized_array.npy"), optimized_array)
        pbar.update(1)

pbar.close()
