# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import numpy as np
import os
import sys
from multiprocessing import Pool, cpu_count, active_children
from SALib.analyze import sobol
from SALib.sample import saltelli
from SALib.analyze import rbd_fast
import pandas as pd
from datetime import datetime
from collections import Counter
# Add project root to Python path so main_code module can be imported
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #


def get_c_well(l_tot, d):
    """

        Return the expected CAPEX for well drilling [€]

         Cost Correlation from:

            Adams et al, “Estimating the Geothermal Electricity Generation Potential of Sedimentary Basins Using
            genGEO (the generalizable GEOthermal techno-economic simulator)", ChemRxiv Prepr., 2021.

    """

    c_well = (0.105 * l_tot ** 2 + 1776 * l_tot * d + 2.753 * pow(10, 5))

    return c_well


def set_parameters_and_execute_avg(k, rho_cp, grad_t, u_inf, mass_flow, depth, l_horiz):
    # Lazy imports: load only when function is called

    # u is the log10 of undisturbed velocity of the flow underground in cm/y
    u = 10**u_inf #from log10(u) to u in cm/y
    u = u/365/24/3600 # from cm/y to cm/s
    u = u/100  # Convert from cm/s to m/s


    from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

        REELWELLHeatingSection,
        REELWELLGeometry,
        REELWEELBHE,

    )
    from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

    days = [0.1, 1, 7, 14, 30, 60, 90, 180, 365, 365*2, 365*5, 365*10, 365*20, 365*30]
    years = (np.array(days)/365).tolist()
    d_well = 0.24

    t_rock = 15 + grad_t * depth
    hs_geometry = REELWELLGeometry(

        l_horiz,
        tub_id=0.085,
        tub_od=0.14,
        cas_id=0.22,
        cas_od=d_well,
        k_insulation=0.1,
        hot_in_tubing=True,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=False,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False

    )
    bhe_in = PlantThermoPoint(["Water"], [1], m_dot=mass_flow)
    bhe_in.set_variable("T", 45)
    bhe_in.set_variable("P", 1.3)

    pe = u * rho_cp * d_well/k
    hs_geometry.rocks_info.pe = pe  # Peclet number now lives in REELWELLRocksInfo, not in REELWEELBHE.__init__

    well = REELWEELBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=15,
        k_rocks=k, c_rocks=1, rho_rocks=rho_cp/1000,
        rw_geometry=hs_geometry, max_iteration=20

    )

    heating_section = REELWELLHeatingSection(

        well, hs_geometry,
        integrate_temperature=True,

    )
    well.heating_section = heating_section
    time = np.array(days)
    w_out_list = []
    t_out_list = []
    p_out_list = []
    time_list = []
    dt_list = []
    t_previous = 0.
    e_tot_avg_kWh = 0
    w_avg = 0
    for t in time:
        well.heating_section.time = t/365
        well.update()
        time_list.append(t)
        t_out_list.append(well.points[-1].get_variable("T"))
        p_out_list.append(well.points[-1].get_variable("P"))
        w_out_list.append(well.power)

        dt_list.append((time_list[-1] - t_previous) * 365 * 24) # in hours

        if t == time[0]:
            if t_out_list[-1] <= 85:
                w_avg = 0
            else:
                w_avg_step = w_out_list[-1]

        else:
            w_avg_step = (w_out_list[-1] + w_previous)/2
        w_previous = w_out_list[-1]
        e_tot_avg_kWh += w_avg_step * dt_list[-1]
        w_avg = e_tot_avg_kWh / (time_list[-1] * 365 * 24)




    if w_avg <= 0.:
        LCOH_ext = 10
        c_well = np.inf
    else:
        c_well = get_c_well(l_horiz+depth, hs_geometry.cas_od) # Add 10% contingency to CAPEX
        LCOH_ext = 1.05 * c_well / (w_avg * time *24 * 0.8)
        i = 0.04
        L_e = years
        CF = 0.8
        alpha = (1 - (1 + i)**(-L_e))/i
        beta = (1+ alpha*1.05)/(alpha * CF * 8760)
        LCOH = beta * c_well / (w_avg)




    return w_avg, LCOH_ext, LCOH


def worker(params_with_idx):
    """
    Worker function to evaluate set_parameters_and_execute for a single parameter combination.
    Uses lazy imports to minimize memory footprint in multiprocessing workers.

    Parameters:
        params_with_idx: tuple (index, params_array)
            - index: int, row index of parameters
            - params_array: array with [k_rock, rho_cp_rock, gradT, pe, mass_flow, depth, l_horiz]

    Returns:
        tuple: (index, w_out, LCOH_ext)
    """
    idx, params = params_with_idx
    try:
        # Ensure project root is in path for imports
        import sys
        import os
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)

        # Lazy imports: only load heavy modules when worker is actually called
        from main_code.well_model.geometry_based_well_models.REELWEEL_model import (
            REELWELLHeatingSection,
            REELWELLGeometry,
            REELWEELBHE,
        )
        from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

        k_rock, rho_cp_rock, gradT, u_inf, mass_flow, depth, l_horiz= params
        w_out, LCOH, lc = set_parameters_and_execute(
            k=k_rock,
            rho_cp=rho_cp_rock,
            grad_t=gradT,
            u_inf=u_inf,
            mass_flow=mass_flow,
            depth=depth,
            l_horiz=l_horiz,
        )
        return idx, w_out, LCOH, lc
    except Exception as e:
        print(f"Worker failed for params index {idx}: {repr(e)}")
        import traceback
        traceback.print_exc()
        # Return default values if computation fails
        return idx, 1e-6, 1e20,1e20


def run_parallel_with_timeout(params_with_idx, n_workers, csv_path, save_every=50,
                              timeout_s=600, timeout_log_path=None):
    """
    Evaluate (index, params) tuples in parallel with a per-point timeout (requires `pebble`).

    A point that exceeds `timeout_s` seconds of execution is killed (its worker process is
    terminated and replaced), saved with the default failure values (1e-6, 1e20, 1e20) and its
    index is appended to `timeout_log_path` so it can be identified/retried later.
    Results are appended to `csv_path` every `save_every` evaluations.
    """
    from pebble import ProcessPool
    from concurrent.futures import as_completed, TimeoutError as FutTimeout
    import gc

    if timeout_log_path is None:
        timeout_log_path = os.path.splitext(csv_path)[0] + '_timeouts.txt'

    results_buffer = []
    n_timeout = 0
    n_error = 0
    pbar = tqdm(total=len(params_with_idx), desc="Evaluating parameter combinations")
    pool = ProcessPool(max_workers=n_workers, max_tasks=100)  # restart workers every 100 tasks
    try:
        futures = {
            pool.schedule(worker, args=(item,), timeout=timeout_s): item[0]
            for item in params_with_idx
        }
        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                result = fut.result()
            except FutTimeout:
                n_timeout += 1
                pbar.write(f"  ⏱ Timeout (>{timeout_s}s) on index {idx}: saved as failed")
                with open(timeout_log_path, 'a') as f:
                    f.write(f"{idx}\n")
                result = (idx, 1e-6, 1e20, 1e20)
            except Exception as e:
                n_error += 1
                pbar.write(f"  ✗ Error on index {idx}: {repr(e)}: saved as failed")
                result = (idx, 1e-6, 1e20, 1e20)

            results_buffer.append(result)
            pbar.update(1)

            if len(results_buffer) >= save_every:
                append_results_to_csv(csv_path, results_buffer)
                results_buffer = []
                gc.collect()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user!")
        raise
    finally:
        if results_buffer:
            append_results_to_csv(csv_path, results_buffer)
            print(f"  ✓ Saved {len(results_buffer)} pending results")
        pbar.close()
        try:
            pool.stop()
            pool.join()
        except Exception:
            pass
        gc.collect()

    print(f"\n  Timeouts: {n_timeout} | Other errors: {n_error}")
    if n_timeout:
        print(f"  Timed-out indices logged in: {timeout_log_path}")


def wrapped_model_parallel(parametres, n_workers=None, csv_path='sobol_results_incremental.csv', save_every=50,
                           timeout_s=600):
    """
    Parallel evaluation of all parameter combinations using multiprocessing with incremental saving.
    WARNING: This may cause memory issues on Windows systems with heavy dependencies.

    Parameters:
        parametres: array of parameter combinations from sobol.sample()
        n_workers: number of parallel workers (default: cpu_count() - 1)
        csv_path: str, path to save incremental results
        save_every: int, save results every N evaluations
        timeout_s: int, max seconds per point; points exceeding it are killed and saved as failed

    Returns:
        tuple: (res_w, res_l) - arrays of results for w_out, and LCOH
    """
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    print(f"Starting parallel evaluation with {n_workers} workers...")
    print(f"Results will be saved every {save_every} iterations to: {csv_path}")

    # Check for existing results to resume from
    completed_indices = load_completed_indices(csv_path)
    if completed_indices:
        print(f"Found {len(completed_indices)} already completed evaluations. Resuming...")
    else:
        # Initialize new results file
        initialize_results_csv(csv_path, include_index=True)

    # Prepare parameter list with indices, skipping completed ones
    params_with_idx = [(idx, p) for idx, p in enumerate(parametres) if idx not in completed_indices]

    if not params_with_idx:
        print("All evaluations already completed! Loading results...")
        df_all = pd.read_csv(csv_path)
        df_all = df_all.sort_values('index').reset_index(drop=True)
        return df_all['w_out'].values, df_all['LCOH'].values

    print(f"Evaluating {len(params_with_idx)} remaining parameter combinations...")

    run_parallel_with_timeout(params_with_idx, n_workers, csv_path,
                              save_every=save_every, timeout_s=timeout_s)

    # Load all results and sort by index
    print(f"\nLoading all results from {csv_path}...")
    df_all = pd.read_csv(csv_path)
    df_all = df_all.sort_values('index').reset_index(drop=True)

    # Extract sorted arrays
    res_w = df_all['w_out'].values

    res_l = df_all['LCOH'].values

    return res_w, res_l


def set_parameters_and_execute(k, rho_cp, grad_t, u_inf, mass_flow, depth, l_horiz, lifetime_years=20):
    # Lazy imports: load only when function is called

    # u is the log10 of undisturbed velocity of the flow underground in cm/y
    u = 10**u_inf #from log10(u) to u in cm/y
    u = u/365/24/3600 # from cm/y to cm/s
    u = u/100  # Convert from cm/s to m/s


    from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

        REELWELLHeatingSection,
        REELWELLGeometry,
        REELWEELBHE,

    )
    from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

    years = lifetime_years  # single evaluation horizon [years]
    d_well = 0.24

    t_rock = 15 + grad_t * depth
    hs_geometry = REELWELLGeometry(

        l_horiz,
        tub_id=0.085,
        tub_od=0.14,
        cas_id=0.22,
        cas_od=d_well,
        k_insulation=0.1,
        hot_in_tubing=True,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=False,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False

    )
    bhe_in = PlantThermoPoint(["Water"], [1], m_dot=mass_flow)
    bhe_in.set_variable("T", 45)
    bhe_in.set_variable("P", 1.3)

    pe = u * rho_cp * d_well/k
    hs_geometry.rocks_info.pe = pe  # Peclet number now lives in REELWELLRocksInfo, not in REELWEELBHE.__init__

    well = REELWEELBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=15,
        k_rocks=k, c_rocks=1, rho_rocks=rho_cp/1000,
        rw_geometry=hs_geometry, max_iteration=20

    )

    heating_section = REELWELLHeatingSection(

        well, hs_geometry,
        integrate_temperature=True,

    )
    well.heating_section = heating_section
    time = 365*years  # [days]

    well.heating_section.time = years  # [years]
    well.update()
    w_out = well.power

    if w_out == 0.:
        LCOH_ext = 10
        c_well = np.inf
        LCOH = 10
    else:
        c_well = get_c_well(l_horiz+depth, hs_geometry.cas_od) # Add 10% contingency to CAPEX

        LCOH_ext = 1.05 * c_well / (w_out * time * 24 * 0.8)
        i = 0.04
        L_e = years
        CF = 0.8
        alpha = (1 - (1 + i)**(-L_e))/i
        beta = (1+ alpha*1.05)/(alpha * CF * 8760)
        LCOH = beta * c_well / (w_out)




    return w_out, LCOH_ext, LCOH


def set_parameters_and_execute_cf(k, rho_cp, grad_t, u_inf, mass_flow, depth, l_horiz, cf):
    # Lazy imports: load only when function is called

    # u is the log10 of undisturbed velocity of the flow underground in cm/y
    u = 10**u_inf #from log10(u) to u in cm/y
    u = u/365/24/3600 # from cm/y to cm/s
    u = u/100  # Convert from cm/s to m/s


    from main_code.well_model.geometry_based_well_models.REELWEEL_model import (

        REELWELLHeatingSection,
        REELWELLGeometry,
        REELWEELBHE,

    )
    from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
    d_well = 0.24
    t_rock = 15 + grad_t * depth
    hs_geometry = REELWELLGeometry(

        l_horiz,
        tub_id=0.085,
        tub_od=0.14,
        cas_id=0.22,
        cas_od=d_well,
        k_insulation=0.1,
        hot_in_tubing=True,
        max_back_time=4,
        alpha_old=0.5,
        neglect_internal_heat_transfer=False,
        ignore_tubing_pressure_losses=False,
        ignore_annulus_pressure_losses=False

    )
    bhe_in = PlantThermoPoint(["Water"], [1], m_dot=mass_flow)
    bhe_in.set_variable("T", 45)
    bhe_in.set_variable("P", 1.3)

    pe = u * rho_cp * d_well/k
    hs_geometry.rocks_info.pe = pe  # Peclet number now lives in REELWELLRocksInfo, not in REELWEELBHE.__init__

    well = REELWEELBHE(

        bhe_in, dz_well=depth, t_rocks=t_rock, t_surf=15,
        k_rocks=k, c_rocks=1, rho_rocks=rho_cp/1000,
        rw_geometry=hs_geometry, max_iteration=20

    )

    heating_section = REELWELLHeatingSection(

        well, hs_geometry,
        integrate_temperature=True,

    )
    well.heating_section = heating_section
    time = 3650*2
    well.heating_section.time = time / 365
    well.update()
    w_out = well.power

    CF = cf


    if w_out == 0.:
        LCOH_ext = 10
        c_well = np.inf
    else:
        c_well = get_c_well(l_horiz+depth, hs_geometry.cas_od) # Add 10% contingency to CAPEX
        LCOH_ext = 1.05 * c_well / (w_out * time *24 * CF)
        i = 0.04
        L_e = 20
        CF = 0.8
        alpha = (1 - (1 + i)**(-L_e))/i
        beta = (1+ alpha*1.05)/(alpha * CF * 8760)
        LCOH = beta * c_well / (w_out)




    return w_out, LCOH_ext, LCOH



# %%------------   VISUALIZATION FUNCTIONS              ---------------------------------------------------> #
def plot_parameters_distribution(df_params, problem, output_folder='Output'):

    os.makedirs(f'{output_folder}/parameter_distributions', exist_ok=True)

    n_params = len(problem['names'])
    n_cols = 3
    n_rows = int(np.ceil(n_params / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
    axes = axes.flatten()

    for idx, param_name in enumerate(problem['names']):
        ax = axes[idx]
        param_values = df_params[param_name].values

        ax.hist(param_values, bins=64, density=False, alpha=0.8, color='skyblue', edgecolor='black')
        ax.set_xlabel(param_name, fontsize=11, fontweight='bold')
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'{param_name} Distribution', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.tight_layout()
    fig.savefig(f'{output_folder}/parameter_distributions/parameters_distribution.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return fig



def plot_sensitivity_indices(problem, sobol_indices, output_name, metric='S1', output_folder='Output'):
    """
    Plot sensitivity indices (S1, ST, etc.) as a bar chart using seaborn.

    Parameters:
        problem: dict with 'names' key containing parameter names
        sobol_indices: DataFrame or dict with sensitivity indices
        output_name: str, name of the output being analyzed
        metric: str, which metric to plot ('S1', 'ST', 'S2', etc.)

    Returns:
        fig: matplotlib figure object
    """
    if isinstance(sobol_indices, pd.DataFrame):
        data = sobol_indices.copy()
    else:
        data = pd.DataFrame(sobol_indices)

    # Extract metric column
    if metric in data.columns:
        values = data[metric].values
    else:
        values = data.iloc[:, 0].values if len(data) > 0 else sobol_indices[metric]

    if "groups" in problem:
        parameters = problem["groups"]
        parameters = list(set(parameters))
    else:
        parameters = problem["names"]
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'Parameter': parameters,
        metric: values
    })

    # Sort by metric value
    plot_data = plot_data.sort_values(metric, ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=plot_data, x=metric, y='Parameter', palette='viridis', ax=ax)
    ax.set_title(f'Sobol {metric} - {output_name}', fontsize=14, fontweight='bold')
    ax.set_xlabel(f'{metric} Index', fontsize=12)
    ax.set_ylabel('Parameters', fontsize=12)
    plt.tight_layout()

    return fig


def plot_sensitivity_heatmap(problem, S1_indices, ST_indices, output_name, output_folder='Output'):
    """
    Plot S1 and ST indices as a heatmap using seaborn.

    Parameters:
        problem: dict with 'names' key containing parameter names
        S1_indices: array-like of first-order sensitivity indices
        ST_indices: array-like of total-order sensitivity indices
        output_name: str, name of the output being analyzed

    Returns:
        fig: matplotlib figure object
    """
    # Create matrix with S1 and ST
    sensitivity_matrix = np.vstack([S1_indices, ST_indices])

    if "groups" in problem:
        parameters = problem["groups"]
        parameters = list(set(parameters))
    else:
        parameters = problem["names"]


    # Create DataFrame for heatmap
    heatmap_data = pd.DataFrame(
        sensitivity_matrix,
        columns=parameters,
        index=['S1 (First-Order)', 'ST (Total-Order)']
    )

    fig, ax = plt.subplots(figsize=(12, 4))
    sns.heatmap(heatmap_data, annot=True, fmt='.4f', cmap='RdYlGn',
                cbar_kws={'label': 'Sensitivity Index'}, ax=ax, vmin=0, vmax=1)
    ax.set_title(f'Sobol Sensitivity Indices - {output_name}', fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def plot_output_2d_map(parametres, problem, output_values, var_idx_1, var_idx_2, output_name, output_folder='Output'):
    """
    Create a 2D heatmap showing output values as a function of two parameters with interpolation.

    Parameters:
        parametres: array of parameter combinations from sobol.sample()
        problem: dict with 'names' key containing parameter names
        output_values: array of output values corresponding to each parameter combination
        var_idx_1: int, index of first parameter to plot
        var_idx_2: int, index of second parameter to plot
        output_name: str, name of the output being analyzed

    Returns:
        fig: matplotlib figure object
    """
    from scipy.interpolate import griddata


    var_name_1 = problem['names'][var_idx_1]
    var_name_2 = problem['names'][var_idx_2]

    # Extract the two parameters and output values
    x = parametres[:, var_idx_1]
    y = parametres[:, var_idx_2]
    z = output_values

    # Remove NaN or infinite values
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x = x[mask]
    y = y[mask]
    z = z[mask]

    if len(x) < 4:
        # Not enough points for interpolation, use simple pivot table
        df = pd.DataFrame(parametres, columns=problem['names'])
        df['output'] = output_values
        pivot_data = df.pivot_table(
            values='output',
            index=var_name_2,
            columns=var_name_1,
            aggfunc='mean'
        )

        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(pivot_data, cmap='viridis', annot=False, cbar_kws={'label': output_name}, ax=ax)
        ax.set_title(f'{output_name} vs {var_name_1} and {var_name_2}',
                     fontsize=14, fontweight='bold')
        ax.set_xlabel(var_name_1, fontsize=12)
        ax.set_ylabel(var_name_2, fontsize=12)
        plt.tight_layout()
        return fig

    # Create a regular grid for interpolation
    grid_resolution = 100  # Increase for smoother interpolation
    xi = np.linspace(x.min(), x.max(), grid_resolution)
    yi = np.linspace(y.min(), y.max(), grid_resolution)
    xi_grid, yi_grid = np.meshgrid(xi, yi)

    # Interpolate using cubic method (smooth), fall back to linear if needed
    try:
        zi_grid = griddata((x, y), z, (xi_grid, yi_grid), method='cubic')
        # Fill any remaining NaN with linear interpolation
        if np.isnan(zi_grid).any():
            zi_linear = griddata((x, y), z, (xi_grid, yi_grid), method='linear')
            zi_grid = np.where(np.isnan(zi_grid), zi_linear, zi_grid)
    except Exception:
        # Fall back to linear interpolation
        zi_grid = griddata((x, y), z, (xi_grid, yi_grid), method='linear')

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Use contourf for filled contour plot (smooth interpolation)
    contour = ax.contourf(xi_grid, yi_grid, zi_grid, levels=20, cmap='viridis')
    cbar = plt.colorbar(contour, ax=ax, label=output_name)

    # Optionally add scatter points to show actual data locations
    scatter = ax.scatter(x, y, c=z, s=10, cmap='viridis', edgecolors='white',
                        linewidth=0.5, alpha=0.6, zorder=5)

    ax.set_title(f'{output_name} vs {var_name_1} and {var_name_2} (Interpolated)',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel(var_name_1, fontsize=12)
    ax.set_ylabel(var_name_2, fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    plt.tight_layout()

    return fig


def plot_all_sensitivity_results(problem, sobol_indices_list, output_names, output_folder='Output'):
    """
    Create a comprehensive comparison plot of all outputs' sensitivity indices.

    Parameters:
        problem: dict with 'names' key containing parameter names
        sobol_indices_list: list of sobol analysis results
        output_names: list of output names (e.g., ['w_out', 'LCOH'])

    Returns:
        fig: matplotlib figure object
    """
    # Collect S1 values for all outputs
    all_S1 = []
    for si in sobol_indices_list:
        # si.to_df() returns (ST_df, S1_df, S2_df), so [1] gets S1
        ST_df, S1_df, S2_df = si.to_df()
        # S1_df is indexed by parameter name, values are in first column or use index
        if 'S1' in S1_df.columns:
            all_S1.append(S1_df['S1'].values)
        else:
            # The values are in the first column, parameter names are in index
            all_S1.append(S1_df.iloc[:, 0].values)

    S1_matrix = np.array(all_S1).T
    if "groups" in problem:
        parameters = problem["groups"]
        parameters = list(set(parameters))
    else:
        parameters = problem["names"]
    # Create DataFrame for heatmap
    heatmap_data = pd.DataFrame(
        S1_matrix,
        columns=output_names,
        index=parameters
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, fmt='.4f', cmap='YlGnBu',
                cbar_kws={'label': 'S1 Index'}, ax=ax, vmin=0, vmax=1)
    ax.set_title('S1 Sensitivity Indices - All Outputs Comparison',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Outputs', fontsize=12)
    ax.set_ylabel('Parameters', fontsize=12)
    plt.tight_layout()

    return fig


# %%------------   SECONDARY PARAMETER GRAPHING FUNCTIONS   ---------------------------------------------------> #

def plot_parameter_correlations(parametres, problem, outputs_dict, output_folder='Output'):
    """
    Create correlation matrix between parameters and outputs.

    Parameters:
        parametres: array of parameter combinations
        problem: dict with 'names' key containing parameter names
        outputs_dict: dict with output names as keys and arrays as values


    Returns:
        fig: matplotlib figure object
    """


    print("\nGenerating parameter correlation matrix...")

    # Create folder for secondary plots
    os.makedirs('Output/secondary_parameter_plots', exist_ok=True)

    try:
        # Combine parameters and outputs into DataFrame
        df_data = pd.DataFrame(parametres, columns=problem['names'])
        for output_name, output_values in outputs_dict.items():
            df_data[output_name] = output_values

        # Calculate correlation matrix
        corr_matrix = df_data.corr()

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 12))
        sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='coolwarm', center=0,
                    square=True, linewidths=0.5, cbar_kws={"shrink": 0.8}, ax=ax)
        ax.set_title('Parameter-Output Correlation Matrix', fontsize=14, fontweight='bold')
        plt.tight_layout()

        # Save figure
        os.makedirs(f'{output_folder}/parameter_correlations', exist_ok=True)
        fig.savefig(f"{output_folder}/parameter_correlations/correlation_matrix.png", dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved correlation matrix")

        return fig

    except Exception as e:
        print(f"  ✗ Failed to generate correlation matrix: {repr(e)}")
        return None


def plot_parameter_vs_output_scatter(parametres, problem, outputs_dict, output_folder='Output'):
    """
    Create scatter plots of each parameter vs each output.

    Parameters:
        parametres: array of parameter combinations
        problem: dict with 'names' key containing parameter names
        outputs_dict: dict with output names as keys and arrays as values


    Returns:
        dict of matplotlib figure objects
    """


    print("\nGenerating parameter vs output scatter plots...")

    # Create folder for scatter plots
    os.makedirs(f'{output_folder}/secondary_parameter_plots/scatter_plots', exist_ok=True)

    figures = {}

    for output_name, output_values in outputs_dict.items():
        # Filter valid values
        valid_mask = np.isfinite(output_values)

        try:
            n_params = len(problem['names'])
            n_cols = 3
            n_rows = int(np.ceil(n_params / n_cols))

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
            axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes

            fig.suptitle(f'Parameter vs {output_name} Scatter Plots', fontsize=16, fontweight='bold')

            for idx, param_name in enumerate(problem['names']):
                ax = axes[idx]
                param_values = parametres[:, idx][valid_mask]
                output_vals = output_values[valid_mask]

                # Create scatter plot with density coloring
                scatter = ax.hexbin(param_values, output_vals, gridsize=30, cmap='viridis', mincnt=1)
                ax.set_xlabel(param_name, fontsize=10)
                ax.set_ylabel(output_name, fontsize=10)
                ax.set_title(param_name, fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3)
                plt.colorbar(scatter, ax=ax, label='Count')

            # Hide unused subplots
            for idx in range(n_params, len(axes)):
                axes[idx].set_visible(False)

            plt.tight_layout()

            # Save figure
            fig.savefig(f"{output_folder}/secondary_parameter_plots/scatter_plots/scatter_{output_name}.png",
                        dpi=300, bbox_inches='tight')
            figures[f"{output_name}_scatter"] = fig
            print(f"  ✓ Saved scatter plots for {output_name}")
            plt.close(fig)

        except Exception as e:
            print(f"  ✗ Failed to generate scatter plots for {output_name}: {repr(e)}")

    return figures


def plot_sensitivity_rankings(problem, sobol_analysis_results, output_folder='Output'):
    """
    Create bar chart showing parameter rankings by sensitivity for all outputs.

    Parameters:
        problem: dict with 'names' key containing parameter names
        sobol_analysis_results: dict from perform_sobol_analysis_from_csv()


    Returns:
        fig: matplotlib figure object
    """


    print("\nGenerating sensitivity rankings plot...")

    # Create folder for secondary plots
    os.makedirs('Output/secondary_parameter_plots', exist_ok=True)

    try:
        sobol_results = sobol_analysis_results['results']

        # Collect ST values (total-order) for all outputs
        all_ST = {}
        for output_name, si in sobol_results.items():
            ST_df, S1_df, S2_df = si.to_df()
            if 'ST' in ST_df.columns:
                all_ST[output_name] = ST_df['ST'].values
            else:
                all_ST[output_name] = ST_df.iloc[:, 0].values

        # Create DataFrame
        df_rankings = pd.DataFrame(all_ST, index=problem['names'])

        # Create figure with grouped bars
        fig, ax = plt.subplots(figsize=(14, 8))

        x = np.arange(len(problem['names']))
        width = 0.8 / len(sobol_results)

        colors = plt.cm.Set3(np.linspace(0, 1, len(sobol_results)))

        for idx, (output_name, values) in enumerate(df_rankings.items()):
            offset = width * idx - (width * (len(sobol_results) - 1) / 2)
            ax.bar(x + offset, values, width, label=output_name, color=colors[idx], alpha=0.8)

        ax.set_xlabel('Parameters', fontsize=12, fontweight='bold')
        ax.set_ylabel('Total Sensitivity Index (ST)', fontsize=12, fontweight='bold')
        ax.set_title('Parameter Sensitivity Rankings (Total-Order Indices)', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(problem['names'], rotation=45, ha='right')
        ax.legend(title='Outputs', loc='upper right')
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()

        # Save figure
        fig.savefig(f"secondary_parameter_plots/sensitivity_rankings.png",
                    dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved sensitivity rankings plot")

        return fig

    except Exception as e:
        print(f"  ✗ Failed to generate sensitivity rankings: {repr(e)}")
        return None


def plot_output_vs_output_scatter(outputs_dict):
    """
    Create scatter plots showing relationships between different outputs.

    Parameters:
        outputs_dict: dict with output names as keys and arrays as values

    Returns:
        fig: matplotlib figure object
    """

    print("\nGenerating output vs output scatter plots...")

    # Create folder for secondary plots
    os.makedirs('Output/secondary_parameter_plots', exist_ok=True)

    try:
        output_names = list(outputs_dict.keys())
        n_outputs = len(output_names)

        if n_outputs < 2:
            print("  ⓘ Need at least 2 outputs for output-vs-output plots")
            return None

        # Create pairwise scatter plots
        fig, axes = plt.subplots(n_outputs, n_outputs, figsize=(4*n_outputs, 4*n_outputs))

        if n_outputs == 2:
            axes = np.array([[axes[0, 0], axes[0, 1]], [axes[1, 0], axes[1, 1]]])

        fig.suptitle('Output Relationships', fontsize=16, fontweight='bold')

        for i, output_name_y in enumerate(output_names):
            for j, output_name_x in enumerate(output_names):
                ax = axes[i, j]

                if i == j:
                    # Diagonal: histogram
                    valid_values = outputs_dict[output_name_x][np.isfinite(outputs_dict[output_name_x])]
                    ax.hist(valid_values, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
                    ax.set_ylabel('Frequency', fontsize=9)
                    if i == n_outputs - 1:
                        ax.set_xlabel(output_name_x, fontsize=10, fontweight='bold')
                else:
                    # Off-diagonal: scatter plot
                    values_x = outputs_dict[output_name_x]
                    values_y = outputs_dict[output_name_y]

                    # Filter valid values
                    valid_mask = np.isfinite(values_x) & np.isfinite(values_y)
                    values_x_valid = values_x[valid_mask]
                    values_y_valid = values_y[valid_mask]

                    if len(values_x_valid) > 0:
                        ax.hexbin(values_x_valid, values_y_valid, gridsize=30, cmap='viridis', mincnt=1)

                    if j == 0:
                        ax.set_ylabel(output_name_y, fontsize=10, fontweight='bold')
                    if i == n_outputs - 1:
                        ax.set_xlabel(output_name_x, fontsize=10, fontweight='bold')

                ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        fig.savefig(f"secondary_parameter_plots/output_relationships.png",
                    dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved output relationships plot")

        return fig

    except Exception as e:
        print(f"  ✗ Failed to generate output relationships: {repr(e)}")
        return None


def plot_parameter_ranges_impact(parametres, problem, outputs_dict):
    """
    Create violin plots showing output distributions for different parameter ranges.

    Parameters:
        parametres: array of parameter combinations
        problem: dict with 'names' key containing parameter names
        outputs_dict: dict with output names as keys and arrays as values

    Returns:
        dict of matplotlib figure objects
    """


    print("\nGenerating parameter ranges impact plots...")

    # Create folder for secondary plots
    os.makedirs('Output/secondary_parameter_plots/range_impact', exist_ok=True)

    figures = {}

    for output_name, output_values in outputs_dict.items():
        try:
            # Filter valid values
            valid_mask = np.isfinite(output_values)

            n_params = len(problem['names'])
            n_cols = 2
            n_rows = int(np.ceil(n_params / n_cols))

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4*n_rows))
            axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes

            fig.suptitle(f'Parameter Ranges Impact on {output_name}', fontsize=16, fontweight='bold')

            for idx, param_name in enumerate(problem['names']):
                ax = axes[idx]
                param_values = parametres[:, idx][valid_mask]
                output_vals = output_values[valid_mask]

                # Split parameter into quartiles
                q1, q2, q3 = np.percentile(param_values, [25, 50, 75])

                # Create bins
                bins = [-np.inf, q1, q2, q3, np.inf]
                labels = ['Low\n(<25%)', 'Med-Low\n(25-50%)', 'Med-High\n(50-75%)', 'High\n(>75%)']

                # Categorize data
                param_categories = pd.cut(param_values, bins=bins, labels=labels)

                # Create DataFrame for violin plot
                df_violin = pd.DataFrame({
                    'Parameter Range': param_categories,
                    output_name: output_vals
                })

                # Violin plot
                sns.violinplot(data=df_violin, x='Parameter Range', y=output_name, ax=ax, palette='Set2')
                ax.set_title(param_name, fontsize=11, fontweight='bold')
                ax.set_xlabel('Parameter Range', fontsize=9)
                ax.set_ylabel(output_name, fontsize=9)
                ax.grid(True, alpha=0.3, axis='y')

            # Hide unused subplots
            for idx in range(n_params, len(axes)):
                axes[idx].set_visible(False)

            plt.tight_layout()

            # Save figure
            fig.savefig(f"secondary_parameter_plots/range_impact/range_impact_{output_name}.png",
                        dpi=300, bbox_inches='tight')
            figures[f"{output_name}_range_impact"] = fig
            print(f"  ✓ Saved parameter range impact for {output_name}")
            plt.close(fig)

        except Exception as e:
            print(f"  ✗ Failed to generate range impact for {output_name}: {repr(e)}")

    return figures


def generate_secondary_plots(problem, parametres, sobol_analysis_results):
    """
    Generate all secondary parameter plots.

    Parameters:
        problem: SALib problem dict
        parametres: array of parameter combinations
        sobol_analysis_results: dict from perform_sobol_analysis_from_csv()

    Returns:
        dict of matplotlib figure objects
    """


    outputs_dict = sobol_analysis_results['outputs']
    figures = {}

    print("\n" + "=" * 80)
    print("GENERATING SECONDARY PARAMETER PLOTS")
    print("=" * 80)

    # 1. Correlation matrix
    fig_corr = plot_parameter_correlations(parametres, problem, outputs_dict)
    if fig_corr is not None:
        figures['correlation_matrix'] = fig_corr
        plt.close(fig_corr)

    # 2. Scatter plots
    scatter_figs = plot_parameter_vs_output_scatter(parametres, problem, outputs_dict)
    figures.update(scatter_figs)

    # 3. Sensitivity rankings
    fig_rankings = plot_sensitivity_rankings(problem, sobol_analysis_results)
    if fig_rankings is not None:
        figures['sensitivity_rankings'] = fig_rankings
        plt.close(fig_rankings)

    # 4. Output relationships
    fig_output_rel = plot_output_vs_output_scatter(outputs_dict)
    if fig_output_rel is not None:
        figures['output_relationships'] = fig_output_rel
        plt.close(fig_output_rel)

    # 5. Parameter ranges impact
    range_figs = plot_parameter_ranges_impact(parametres, problem, outputs_dict)
    figures.update(range_figs)

    print(f"\n✓ Generated {len(figures)} secondary parameter plots")
    return figures


# %%------------   SOBOL ANALYSIS FUNCTIONS (Standalone)      ---------------------------------------------------> #

def perform_sobol_analysis_from_csv(
    csv_results_file,
    csv_params_file,
    problem,
    output_columns=None,
    lcoh_max=None,
    w_min=None
):
    """
    Perform Sobol sensitivity analysis from CSV results, checking completed indices.

    Parameters:
        csv_results_file: path to results CSV (expects 'index' column when available)
        csv_params_file: path to parameters CSV
        problem: SALib problem dict
        output_columns: list of outputs to analyze (default: ['w_out', 'LCOH'])
        lcoh_max: optional maximum LCOH threshold (values above this are excluded)
        w_min: optional minimum w_effective threshold (values below this are excluded)

    Returns:
        dict with analysis results and completion metadata
    """
    if output_columns is None:
        output_columns = ['w_out', 'LCOH']


    print("\n" + "=" * 80)
    print("SOBOL SENSITIVITY ANALYSIS FROM CSV RESULTS")
    print("=" * 80)

    print(f"\nLoading parameters from: {csv_params_file}")
    df_params = pd.read_csv(csv_params_file)
    parametres = df_params.values
    print(f"  ✓ Loaded {len(parametres)} total parameter combinations")

    print(f"\nLoading results from: {csv_results_file}")
    df_results = pd.read_csv(csv_results_file)
    print(f"  ✓ Loaded {len(df_results)} results")

    # Check for 'index' column and extract completed indices
    if 'index' not in df_results.columns:
        print("  ⚠ 'index' column not found in results; using row order.")
        # Assume results are in order 0, 1, 2, ... N-1
        completed_indices = set(range(len(df_results)))
        completed_indices_sorted = list(range(len(df_results)))
    else:
        # Extract indices from results, sort, and filter parameters accordingly
        completed_indices = set(df_results['index'].astype(int).values)
        df_results = df_results.sort_values('index').reset_index(drop=True)
        completed_indices_sorted = df_results['index'].astype(int).values.tolist()

    all_indices = set(range(len(parametres)))
    missing_indices = sorted(list(all_indices - completed_indices))

    print(f"  ✓ Completed evaluations: {len(completed_indices)}/{len(parametres)}")
    if missing_indices:
        print(f"  ⚠ Missing {len(missing_indices)} evaluations")
        if len(missing_indices) <= 20:
            print(f"     Missing indices: {missing_indices}")
        else:
            print(f"     Missing indices (first 20): {missing_indices[:20]}")

    # Filter parameters to match only the completed indices
    parametres_completed = parametres[completed_indices_sorted]

    # Verify alignment
    if len(parametres_completed) != len(df_results):
        print(f"  ✗ ERROR: Parameter count ({len(parametres_completed)}) != Result count ({len(df_results)})")
        return None

    print(f"\n  ✓ Using {len(parametres_completed)} completed evaluations for analysis")
    print(f"  ✓ Parameters and results are properly aligned")

    # Extract output columns from results
    outputs_dict = {}
    for col in output_columns:
        if col in df_results.columns:
            outputs_dict[col] = df_results[col].values
            print(f"  ✓ Extracted '{col}' output ({len(df_results[col].values)} values)")
        else:
            print(f"  ✗ Column '{col}' not found in results")

    if not outputs_dict:
        print("  ✗ No valid output columns found; aborting analysis")
        return None

    # Apply value filters if specified
    print("\n" + "=" * 80)
    print("APPLYING OUTPUT VALUE FILTERS")
    print("=" * 80)

    filter_mask = np.ones(len(df_results), dtype=bool)
    filtered_count = 0

    if lcoh_max is not None and 'LCOH' in outputs_dict:
        lcoh_filter = outputs_dict['LCOH'] <= lcoh_max
        excluded_lcoh = np.sum(~lcoh_filter)
        if excluded_lcoh > 0:
            print(f"\nLCOH Filter: max = {lcoh_max:.6e}")
            print(f"  ✗ Excluded {excluded_lcoh} values above threshold")
            filter_mask &= lcoh_filter
            filtered_count += excluded_lcoh
    else:
        print(f"\nLCOH Filter: None (keeping all values)")

    if w_min is not None and 'w_out' in outputs_dict:
        w_filter = outputs_dict['w_out'] >= w_min
        excluded_w = np.sum(~w_filter)
        if excluded_w > 0:
            print(f"\nw_out Filter: min = {w_min:.6e}")
            print(f"  ✗ Excluded {excluded_w} values below threshold")
            filter_mask &= w_filter
            filtered_count += excluded_w
    else:
        print(f"\nw_out Filter: None (keeping all values)")

    # Apply filter to all outputs and parameters
    if filtered_count > 0:
        print(f"\n  Total rows excluded by filters: {filtered_count}")
        print(f"  Remaining rows for analysis: {np.sum(filter_mask)}")

        for col in list(outputs_dict.keys()):
            outputs_dict[col] = outputs_dict[col][filter_mask]

        parametres_completed = parametres_completed[filter_mask]
        df_results = df_results[filter_mask].reset_index(drop=True)
    else:
        print(f"\n  ✓ No values filtered; all data retained")

    print(f"\n  ✓ Using {len(parametres_completed)} evaluations for analysis")
    print("\n" + "=" * 80)
    print("OUTPUT VALUE STATISTICS")
    print("=" * 80)
    for col, values in outputs_dict.items():
        valid_values = values[np.isfinite(values)]
        if len(valid_values) > 0:
            print(f"\n{col}:")
            print(f"  Count: {len(valid_values)}")
            print(f"  Min: {np.min(valid_values):.6e}")
            print(f"  Max: {np.max(valid_values):.6e}")
            print(f"  Mean: {np.mean(valid_values):.6e}")
            print(f"  Median: {np.median(valid_values):.6e}")
            print(f"  Std Dev: {np.std(valid_values):.6e}")
            print(f"  Q1 (25%): {np.percentile(valid_values, 25):.6e}")
            print(f"  Q3 (75%): {np.percentile(valid_values, 75):.6e}")

            # Count outliers (values > 3 std dev from mean)
            mean = np.mean(valid_values)
            std = np.std(valid_values)
            outliers = np.sum(np.abs(valid_values - mean) > 3 * std)
            if outliers > 0:
                print(f"  Outliers (>3σ): {outliers} ({100*outliers/len(valid_values):.2f}%)")

    print("\n" + "=" * 80)
    print("PERFORMING SOBOL ANALYSIS")
    print("=" * 80)

    # Verify SALib requirements
    # SALib needs N*(2D+2) samples where N is base sample size, D is number of parameters
    if "groups" in problem:
        cc = Counter(problem['groups'])
        pattern_size = 2 * len(cc) + 2
    else:
        pattern_size = 2 * problem['num_vars'] + 2
    expected_base_size = len(parametres_completed) / pattern_size

    print(f"\n  Sobol sample info:")
    print(f"    - Total samples available: {len(parametres_completed)}")
    print(f"    - Number of parameters (D): {problem['num_vars']}")
    print(f"    - Expected pattern: N*(2D+2) = N*{pattern_size}")
    print(f"    - Calculated base size (N): {expected_base_size:.2f}")

    # Handle incomplete Sobol sample sets
    if not expected_base_size.is_integer():
        # Calculate nearest complete sample set
        n_complete = int(expected_base_size)  # Round down to nearest integer
        n_samples_needed = n_complete * pattern_size

        print(f"\n  ⚠ WARNING: Sample count doesn't match Sobol pattern!")
        print(f"     Have: {len(parametres_completed)} samples")
        print(f"     Need: Multiple of {pattern_size} (e.g., {(n_complete+1)*pattern_size})")
        print(f"     → Truncating to {n_samples_needed} samples (N={n_complete})")

        # Truncate to complete set
        parametres_completed = parametres_completed[:n_samples_needed]

        # Also truncate outputs
        for col in list(outputs_dict.keys()):
            outputs_dict[col] = outputs_dict[col][:n_samples_needed]
            print(f"     → Truncated '{col}' to {len(outputs_dict[col])} values")

    sobol_results = {}
    for output_name, output_values in outputs_dict.items():
        print(f"\nAnalyzing {output_name} sensitivity...")

        # Double-check alignment
        if len(output_values) != len(parametres_completed):
            print(f"  ✗ ERROR: Output length ({len(output_values)}) != Parameter length ({len(parametres_completed)})")
            continue

        try:
            # Perform Sobol analysis using only the completed parameter-output pairs
            si = sobol.analyze(problem, output_values, print_to_console=True)
            sobol_results[output_name] = si

            si_T, si_1, si_2 = si.to_df()
            si_T.to_csv(f"sobol_analysis_{output_name}_ST.csv", index=True)
            si_1.to_csv(f"sobol_analysis_{output_name}_S1.csv", index=True)
            si_2.to_csv(f"sobol_analysis_{output_name}_S2.csv", index=True)
            print(f"  ✓ Saved analysis for {output_name}")
        except Exception as e:
            print(f"  ✗ Failed to analyze {output_name}: {repr(e)}")
            import traceback
            traceback.print_exc()

    return {
        'results': sobol_results,
        'outputs': outputs_dict,
        'parameters': parametres_completed,
        'completed_count': len(completed_indices),
        'total_count': len(parametres),
        'missing_indices': missing_indices,
    }


def perform_sobol_visualization(
    problem,
    sobol_analysis_results,
    output_columns=None,
    generate_2d_maps=True,
    parametres=None
):
    """
    Generate visualization plots from Sobol analysis results, organizing them into folders.

    Parameters:
        problem: SALib problem dict
        sobol_analysis_results: dict from perform_sobol_analysis_from_csv()
        output_columns: list of outputs to visualize
        generate_2d_maps: whether to generate 2D output maps
        parametres: parameter array (required if generate_2d_maps=True)

    Returns:
        dict of matplotlib figure objects
    """
    if output_columns is None:
        output_columns = ['w_out', 'LCOH']


    if sobol_analysis_results is None:
        print("✗ No analysis results provided")
        return {}

    # Create folders for organizing plots
    folders = {
        'sensitivity_indices': 'sensitivity_indices_barplots',
        'heatmaps': 'sensitivity_heatmaps',
        'comparison': 'sensitivity_comparison',
        'maps_2d': 'output_2d_maps'
    }

    for folder in folders.values():
        os.makedirs(folder, exist_ok=True)

    sobol_results = sobol_analysis_results['results']
    outputs_dict = sobol_analysis_results['outputs']
    figures = {}

    print("\n" + "=" * 80)
    print("GENERATING VISUALIZATION PLOTS")
    print("=" * 80)

    # ...existing code for S1 and ST plots...
    for output_name in output_columns:
        if output_name not in sobol_results:
            print(f"  ⊘ Skipping {output_name} (no analysis results)")
            continue

        si = sobol_results[output_name]
        si_T, si_1, _ = si.to_df()

        try:
            fig_s1 = plot_sensitivity_indices(problem, si_1, output_name, metric='S1')
            fig_s1.savefig(f"{folders['sensitivity_indices']}/S1_barplot_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_S1"] = fig_s1
            plt.close(fig_s1)
        except Exception as e:
            print(f"  ✗ Failed S1 plot for {output_name}: {repr(e)}")

        try:
            fig_st = plot_sensitivity_indices(problem, si_T, output_name, metric='ST')
            fig_st.savefig(f"{folders['sensitivity_indices']}/ST_barplot_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_ST"] = fig_st
            plt.close(fig_st)
        except Exception as e:
            print(f"  ✗ Failed ST plot for {output_name}: {repr(e)}")

    # ...existing code for heatmaps...
    for output_name in output_columns:
        if output_name not in sobol_results:
            continue

        si = sobol_results[output_name]
        si_T, si_1, _ = si.to_df()

        si_1_values = si_1['S1'].values if 'S1' in si_1.columns else si_1.iloc[:, 0].values
        si_T_values = si_T['ST'].values if 'ST' in si_T.columns else si_T.iloc[:, 0].values

        try:
            fig_heat = plot_sensitivity_heatmap(problem, si_1_values, si_T_values, output_name)
            fig_heat.savefig(f"{folders['heatmaps']}/heatmap_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_heatmap"] = fig_heat
            plt.close(fig_heat)
        except Exception as e:
            print(f"  ✗ Failed heatmap for {output_name}: {repr(e)}")

    # ...existing code for comparison...
    if len(sobol_results) > 1:
        try:
            sobol_list = [sobol_results[col] for col in output_columns if col in sobol_results]
            output_names_list = [col for col in output_columns if col in sobol_results]
            fig_all = plot_all_sensitivity_results(problem, sobol_list, output_names_list)
            fig_all.savefig(f"{folders['comparison']}/all_outputs_comparison.png", dpi=300, bbox_inches='tight')
            figures['all_outputs_comparison'] = fig_all
            plt.close(fig_all)
        except Exception as e:
            print(f"  ✗ Failed comparison heatmap: {repr(e)}")

    # ...existing code for 2D maps...
    if generate_2d_maps and parametres is not None:
        param_pairs = [
            (0, 1),
            (0, 4),
            (1, 4),
            (2, 4),
            (4, 5),
            (5, 6),
        ]
        for idx1, idx2 in param_pairs:
            var1 = problem['names'][idx1]
            var2 = problem['names'][idx2]
            for output_name, output_values in outputs_dict.items():
                try:
                    fig_map = plot_output_2d_map(parametres, problem, output_values, idx1, idx2, output_name)
                    fig_map.savefig(
                        f"{folders['maps_2d']}/{output_name}_{var1}_vs_{var2}.png",
                        dpi=300,
                        bbox_inches='tight'
                    )
                    figures[f"{output_name}_map_{var1}_{var2}"] = fig_map
                    plt.close(fig_map)
                except Exception as e:
                    print(f"  ✗ Failed 2D map for {output_name} ({var1} vs {var2}): {repr(e)}")

    print(f"\n✓ Visualization complete: {len(figures)} plots generated")
    return figures


# %%------------   OUTPUT DISTRIBUTION PLOTTING FUNCTIONS              ---------------------------------------------------> #
def plot_output_distributions(outputs_dict):
    """
    Create distribution plots (histogram + KDE) for all outputs with outlier filtering.
    Saves each plot type separately in dedicated folders.

    Parameters:
        outputs_dict: dict with output names as keys and arrays as values

    Returns:
        dict of matplotlib figure objects
    """


    figures = {}

    print("\n" + "=" * 80)
    print("GENERATING OUTPUT DISTRIBUTION PLOTS")
    print("=" * 80)

    # Create folders for each plot type
    folders = {
        'histogram': 'distribution_histograms',
        'boxplot': 'distribution_boxplots',
    }

    for folder in folders.values():
        os.makedirs(folder, exist_ok=True)

    for output_name, output_values in outputs_dict.items():
        print(f"\nGenerating distribution plots for '{output_name}'...")

        try:
            # Remove invalid values for statistics
            valid_values = output_values[np.isfinite(output_values)]

            if len(valid_values) == 0:
                print(f"  ✗ No valid values for {output_name}")
                continue

            # Filter extreme outliers (values > 3 std dev from mean) for better visualization
            mean = np.mean(valid_values)
            std = np.std(valid_values)
            filtered_values = valid_values[np.abs(valid_values - mean) <= 3 * std]

            outlier_count = len(valid_values) - len(filtered_values)
            if outlier_count > 0:
                print(f"  ⓘ Filtered {outlier_count} extreme outliers (>3σ) from {output_name} for visualization")

            # 1. Histogram with KDE
            fig1, ax1 = plt.subplots(figsize=(12, 6))
            ax1.hist(filtered_values, bins=50, density=True, alpha=0.6, color='skyblue', edgecolor='black')
            filtered_values_series = pd.Series(filtered_values)
            filtered_values_series.plot(kind='kde', ax=ax1, color='darkblue', linewidth=2)
            ax1.set_xlabel(output_name, fontsize=11)
            ax1.set_ylabel('Density', fontsize=11)
            ax1.set_title(f'Histogram with KDE - {output_name} (Filtered, n={len(filtered_values)})', fontsize=12, fontweight='bold')
            ax1.grid(True, alpha=0.3)
            plt.tight_layout()
            fig1.savefig(f"{folders['histogram']}/histogram_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_histogram"] = fig1
            print(f"  ✓ Saved histogram for {output_name}")
            plt.close(fig1)

            # 2. Box plot
            fig2, ax2 = plt.subplots(figsize=(10, 6))
            bp = ax2.boxplot(filtered_values, vert=True, patch_artist=True)
            bp['boxes'][0].set_facecolor('lightblue')
            ax2.set_ylabel(output_name, fontsize=11)
            ax2.set_title(f'Box Plot - {output_name} (Filtered)', fontsize=12, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='y')

            # Add statistics text (all values)
            stats_text = (
                f"Count: {len(valid_values)}\n"
                f"Mean: {np.mean(valid_values):.4e}\n"
                f"Median: {np.median(valid_values):.4e}\n"
                f"Std Dev: {np.std(valid_values):.4e}\n"
                f"Min: {np.min(valid_values):.4e}\n"
                f"Max: {np.max(valid_values):.4e}\n"
                f"Q1: {np.percentile(valid_values, 25):.4e}\n"
                f"Q3: {np.percentile(valid_values, 75):.4e}"
            )
            ax2.text(1.3, np.median(filtered_values), stats_text, fontsize=9,
                    verticalalignment='center', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            plt.tight_layout()
            fig2.savefig(f"{folders['boxplot']}/boxplot_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_boxplot"] = fig2
            print(f"  ✓ Saved box plot for {output_name}")
            plt.close(fig2)





        except Exception as e:
            print(f"  ✗ Failed to generate distributions for {output_name}: {repr(e)}")

    return figures


def plot_output_comparison(outputs_dict):
    """
    Create a comparison plot showing all outputs on the same figure with outlier filtering.

    Parameters:
        outputs_dict: dict with output names as keys and arrays as values

    Returns:
        matplotlib figure object
    """


    print(f"\nGenerating comparison plot for all outputs...")

    # Create folder for comparison plots
    os.makedirs('Output/distribution_comparison', exist_ok=True)

    try:
        n_outputs = len(outputs_dict)
        fig, axes = plt.subplots(1, n_outputs, figsize=(6*n_outputs, 5))

        # Handle single output case
        if n_outputs == 1:
            axes = [axes]

        fig.suptitle('Output Distributions Comparison (Filtered)', fontsize=16, fontweight='bold')

        for ax, (output_name, output_values) in zip(axes, outputs_dict.items()):
            valid_values = output_values[np.isfinite(output_values)]

            if len(valid_values) > 0:
                # Filter extreme outliers for visualization
                mean = np.mean(valid_values)
                std = np.std(valid_values)
                filtered_values = valid_values[np.abs(valid_values - mean) <= 3 * std]

                ax.hist(filtered_values, bins=40, density=True, alpha=0.6, color='skyblue', edgecolor='black')
                filtered_values_series = pd.Series(filtered_values)
                filtered_values_series.plot(kind='kde', ax=ax, color='darkblue', linewidth=2)
                ax.set_xlabel(output_name, fontsize=11)
                ax.set_ylabel('Density', fontsize=11)
                ax.set_title(f'{output_name}\n(n={len(filtered_values)})', fontsize=12, fontweight='bold')
                ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        fig.savefig(f"distribution_comparison/comparison.png", dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved comparison plot")

        return fig

    except Exception as e:
        print(f"  ✗ Failed to generate comparison plot: {repr(e)}")
        return None



# %%------------   RBD-FAST ANALYSIS FUNCTIONS

def perform_rbd_fast_analysis_from_csv(
    csv_results_file,
    csv_params_file,
    problem,
    output_columns=None,
    lcoh_max=None,
    w_min=None
):
    """
    Perform RBD-FAST sensitivity analysis from CSV results.

    RBD-FAST (Radial Basis Design - Factorial At Space-filling Technique)
    is an alternative to Sobol for global sensitivity analysis. It requires
    fewer model evaluations but provides similar insights.

    Parameters:
        csv_results_file: path to results CSV (expects 'index' column when available)
        csv_params_file: path to parameters CSV
        problem: SALib problem dict
        output_columns: list of outputs to analyze (default: ['w_out', 'LCOH'])
        lcoh_max: optional maximum LCOH threshold (values above this are excluded)
        w_min: optional minimum w_effective threshold (values below this are excluded)

    Returns:
        dict with analysis results and completion metadata
    """
    if output_columns is None:
        output_columns = ['w_out', 'LCOH', 'LCOH_i']


    print("\n" + "=" * 80)
    print("RBD-FAST SENSITIVITY ANALYSIS FROM CSV RESULTS")
    print("=" * 80)

    print(f"\nLoading parameters from: {csv_params_file}")
    df_params = pd.read_csv(csv_params_file)
    parametres = df_params.values
    print(f"  ✓ Loaded {len(parametres)} total parameter combinations")

    print(f"\nLoading results from: {csv_results_file}")
    df_results = pd.read_csv(csv_results_file)
    print(f"  ✓ Loaded {len(df_results)} results")

    # Check for 'index' column and extract completed indices
    if 'index' not in df_results.columns:
        print("  ⚠ 'index' column not found in results; using row order.")
        completed_indices = set(range(len(df_results)))
        completed_indices_sorted = list(range(len(df_results)))
    else:
        # Extract indices from results, sort, and filter parameters accordingly
        completed_indices = set(df_results['index'].astype(int).values)
        df_results = df_results.sort_values('index').reset_index(drop=True)
        completed_indices_sorted = df_results['index'].astype(int).values.tolist()

    all_indices = set(range(len(parametres)))
    missing_indices = sorted(list(all_indices - completed_indices))

    print(f"  ✓ Completed evaluations: {len(completed_indices)}/{len(parametres)}")
    if missing_indices:
        print(f"  ⚠ Missing {len(missing_indices)} evaluations")
        if len(missing_indices) <= 20:
            print(f"     Missing indices: {missing_indices}")
        else:
            print(f"     Missing indices (first 20): {missing_indices[:20]}")

    # Filter parameters to match only the completed indices
    parametres_completed = parametres[completed_indices_sorted]

    # Verify alignment
    if len(parametres_completed) != len(df_results):
        print(f"  ✗ ERROR: Parameter count ({len(parametres_completed)}) != Result count ({len(df_results)})")
        return None

    print(f"\n  ✓ Using {len(parametres_completed)} completed evaluations for analysis")
    print(f"  ✓ Parameters and results are properly aligned")

    # Extract output columns from results
    outputs_dict = {}
    for col in output_columns:
        if col in df_results.columns:
            outputs_dict[col] = df_results[col].values
            print(f"  ✓ Extracted '{col}' output ({len(df_results[col].values)} values)")
        else:
            print(f"  ✗ Column '{col}' not found in results")

    if not outputs_dict:
        print("  ✗ No valid output columns found; aborting analysis")
        return None

    # Apply value filters if specified
    print("\n" + "=" * 80)
    print("APPLYING OUTPUT VALUE FILTERS")
    print("=" * 80)

    filter_mask = np.ones(len(df_results), dtype=bool)
    filtered_count = 0

    if lcoh_max is not None and 'LCOH' in outputs_dict:
        lcoh_filter = outputs_dict['LCOH'] >= lcoh_max
        excluded_lcoh = np.sum(~lcoh_filter)
        if excluded_lcoh > 0:
            print(f"\nLCOH Filter: max = {lcoh_max:.6e}")
            print(f"  ✗ Excluded {excluded_lcoh} values above threshold")
            filter_mask &= lcoh_filter
            filtered_count += excluded_lcoh
    else:
        print(f"\nLCOH Filter: None (keeping all values)")

    if w_min is not None and 'w_out' in outputs_dict:
        w_filter = outputs_dict['w_out'] >= w_min
        excluded_w = np.sum(~w_filter)
        if excluded_w > 0:
            print(f"\nw_out Filter: min = {w_min:.6e}")
            print(f"  ✗ Excluded {excluded_w} values below threshold")
            filter_mask &= w_filter
            filtered_count += excluded_w
    else:
        print(f"\nw_out Filter: None (keeping all values)")

    # Apply filter to all outputs and parameters
    if filtered_count > 0:
        print(f"\n  Total rows excluded by filters: {filtered_count}")
        print(f"  Remaining rows for analysis: {np.sum(filter_mask)}")

        for col in list(outputs_dict.keys()):
            outputs_dict[col] = outputs_dict[col][filter_mask]

        parametres_completed = parametres_completed[filter_mask]
        df_results = df_results[filter_mask].reset_index(drop=True)
    else:
        print(f"\n  ✓ No values filtered; all data retained")

    print(f"\n  ✓ Using {len(parametres_completed)} evaluations for analysis")
    print("\n" + "=" * 80)
    print("OUTPUT VALUE STATISTICS")
    print("=" * 80)
    for col, values in outputs_dict.items():
        valid_values = values[np.isfinite(values)]
        if len(valid_values) > 0:
            print(f"\n{col}:")
            print(f"  Count: {len(valid_values)}")
            print(f"  Min: {np.min(valid_values):.6e}")
            print(f"  Max: {np.max(valid_values):.6e}")
            print(f"  Mean: {np.mean(valid_values):.6e}")
            print(f"  Median: {np.median(valid_values):.6e}")
            print(f"  Std Dev: {np.std(valid_values):.6e}")
            print(f"  Q1 (25%): {np.percentile(valid_values, 25):.6e}")
            print(f"  Q3 (75%): {np.percentile(valid_values, 75):.6e}")

            # Count outliers (values > 3 std dev from mean)
            mean = np.mean(valid_values)
            std = np.std(valid_values)
            outliers = np.sum(np.abs(valid_values - mean) > 3 * std)
            if outliers > 0:
                print(f"  Outliers (>3σ): {outliers} ({100*outliers/len(valid_values):.2f}%)")

    print("\n" + "=" * 80)
    print("PERFORMING RBD-FAST ANALYSIS")
    print("=" * 80)

    rbd_fast_results = {}
    for output_name, output_values in outputs_dict.items():
        print(f"\nAnalyzing {output_name} sensitivity with RBD-FAST...")

        # Double-check alignment
        if len(output_values) != len(parametres_completed):
            print(f"  ✗ ERROR: Output length ({len(output_values)}) != Parameter length ({len(parametres_completed)})")
            continue

        try:
            # Perform RBD-FAST analysis
            si = rbd_fast.analyze(problem, parametres_completed, output_values, print_to_console=True)
            rbd_fast_results[output_name] = si

            # Save results to CSV
            si_df = si.to_df()
            si_df.to_csv(f"rbd_fast_analysis_{output_name}.csv", index=True)
            print(f"  ✓ Saved RBD-FAST analysis for {output_name}")
        except Exception as e:
            print(f"  ✗ Failed to analyze {output_name} with RBD-FAST: {repr(e)}")
            import traceback
            traceback.print_exc()

    return {
        'results': rbd_fast_results,
        'outputs': outputs_dict,
        'parameters': parametres_completed,
        'completed_count': len(completed_indices),
        'total_count': len(parametres),
        'missing_indices': missing_indices,
    }

def perform_rbd_from_dicts(
    problem,
    parameters_dict,
    outputs_dict,
    string = None
):
    """
    Perform RBD-FAST analysis directly from parameter and output dictionaries.

    Parameters:
        problem: SALib problem dict
        parametres_dict: dict with parameter names as keys and arrays as values
        outputs_dict: dict with output names as keys and arrays as values
    Returns:
        dict with RBD-FAST analysis results
    """
    print("\n" + "=" * 80)
    print("PERFORMING RBD-FAST ANALYSIS FROM DICTIONARIES")
    print("=" * 80)

    rbd_fast_results = {}
    for key in outputs_dict.keys():
        print(f"\nAnalyzing {key} sensitivity with RBD-FAST...")
        csv_name = string + key if string else key
        try:
            # Perform RBD-FAST analysis
            print(key)
            si = rbd_fast.analyze(problem, parameters_dict, outputs_dict[key], print_to_console=True)
            rbd_fast_results[key] = si

            # Save results to CSV
            si_df = si.to_df()
            si_df.to_csv(f"rbd_fast_analysis_{csv_name}.csv", index=True)
            print(f"  ✓ Saved RBD-FAST analysis for {key}")
        except Exception as e:
            print(f"  ✗ Failed to analyze {key} with RBD-FAST: {repr(e)}")
            import traceback
            traceback.print_exc()

    return {
        'results': rbd_fast_results,
        'outputs': outputs_dict,
    }

def perform_rbd_fast_visualization(
    problem,
    rbd_fast_analysis_results,
    output_columns=None,
    parametres=None
):
    """
    Generate visualization plots from RBD-FAST analysis results.

    Parameters:
        problem: SALib problem dict
        rbd_fast_analysis_results: dict from perform_rbd_fast_analysis_from_csv()
        output_columns: list of outputs to visualize
        parametres: parameter array (for 2D maps)

    Returns:
        dict of matplotlib figure objects
    """
    if output_columns is None:
        output_columns = ['w_out', 'LCOH', 'LCOH_i']


    if rbd_fast_analysis_results is None:
        print("✗ No analysis results provided")
        return {}

    # Create folders for organizing plots
    folders = {
        'rbd_indices': 'rbd_fast_barplots',
        'heatmaps': 'rbd_fast_heatmaps',
        'comparison': 'rbd_fast_comparison'
    }

    for folder in folders.values():
        os.makedirs(folder, exist_ok=True)

    rbd_results = rbd_fast_analysis_results['results']
    outputs_dict = rbd_fast_analysis_results['outputs']
    figures = {}

    print("\n" + "=" * 80)
    print("GENERATING RBD-FAST VISUALIZATION PLOTS")
    print("=" * 80)

    # Generate bar plots for each output
    for output_name in output_columns:
        if output_name not in rbd_results:
            print(f"  ⊘ Skipping {output_name} (no analysis results)")
            continue

        si = rbd_results[output_name]
        si_df = si.to_df()

        try:
            # Extract sensitivity index column (usually named 'S1')
            if 'S1' in si_df.columns:
                si_values = si_df['S1'].values
            else:
                si_values = si_df.iloc[:, 0].values

            # Create bar plot
            param_names = problem['names']
            plot_data = pd.DataFrame({
                'Parameter': param_names,
                'S1': si_values
            })
            plot_data = plot_data.sort_values('S1', ascending=True)

            fig, ax = plt.subplots(figsize=(10, 6))
            sns.barplot(data=plot_data, x='S1', y='Parameter', palette='viridis', ax=ax)
            ax.set_title(f'RBD-FAST Sensitivity Index - {output_name}', fontsize=14, fontweight='bold')
            ax.set_xlabel('Sensitivity Index', fontsize=12)
            ax.set_ylabel('Parameters', fontsize=12)
            plt.tight_layout()

            fig.savefig(f"{folders['rbd_indices']}/rbd_fast_barplot_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_rbd_fast"] = fig
            print(f"  ✓ Generated RBD-FAST bar plot for {output_name}")
            plt.close(fig)

        except Exception as e:
            print(f"  ✗ Failed to create bar plot for {output_name}: {repr(e)}")

    # Generate comparison heatmap if multiple outputs
    if len(rbd_results) > 1:
        print("\nGenerating RBD-FAST comparison heatmap...")
        try:
            # Collect all sensitivity indices
            all_si = {}
            for output_name in output_columns:
                if output_name in rbd_results:
                    si = rbd_results[output_name]
                    si_df = si.to_df()
                    if 'S1' in si_df.columns:
                        all_si[output_name] = si_df['S1'].values
                    else:
                        all_si[output_name] = si_df.iloc[:, 0].values

            # Create heatmap
            heatmap_data = pd.DataFrame(all_si, index=problem['names'])

            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(heatmap_data, annot=True, fmt='.4f', cmap='YlGnBu',
                        cbar_kws={'label': 'Sensitivity Index'}, ax=ax)
            ax.set_title('RBD-FAST Sensitivity Indices - All Outputs Comparison',
                         fontsize=14, fontweight='bold')
            ax.set_xlabel('Outputs', fontsize=12)
            ax.set_ylabel('Parameters', fontsize=12)
            plt.tight_layout()

            fig.savefig(f"{folders['comparison']}/rbd_fast_comparison.png", dpi=300, bbox_inches='tight')
            figures['rbd_fast_comparison'] = fig
            print(f"  ✓ Generated RBD-FAST comparison heatmap")
            plt.close(fig)

        except Exception as e:
            print(f"  ✗ Failed to create comparison heatmap: {repr(e)}")

    print(f"\n✓ RBD-FAST Visualization complete: {len(figures)} plots generated")
    return figures


def compare_rbd_fast_vs_sobol(
    sobol_results,
    rbd_fast_results,
    problem,
    output_columns=None,
):
    """
    Create comparison plots between Sobol and RBD-FAST sensitivity analyses.

    Parameters:
        sobol_results: dict from perform_sobol_analysis_from_csv()
        rbd_fast_results: dict from perform_rbd_fast_analysis_from_csv()
        problem: SALib problem dict
        output_columns: list of outputs to compare

    Returns:
        dict of matplotlib figure objects
    """
    if output_columns is None:
        output_columns = ['w_out', 'LCOH']


    print("\n" + "=" * 80)
    print("COMPARING SOBOL VS RBD-FAST SENSITIVITY ANALYSIS")
    print("=" * 80)

    # Create folder for comparison plots
    os.makedirs('Output/method_comparison', exist_ok=True)
    figures = {}

    sobol_si = sobol_results['results']
    rbd_fast_si = rbd_fast_results['results']

    for output_name in output_columns:
        if output_name not in sobol_si or output_name not in rbd_fast_si:
            print(f"  ⊘ Skipping {output_name} (results not available for both methods)")
            continue

        try:
            # Extract sensitivity indices from both methods
            sobol_result = sobol_si[output_name]
            rbd_result = rbd_fast_si[output_name]

            # Sobol: need to extract S1 values
            sobol_ST, sobol_S1, sobol_S2 = sobol_result.to_df()
            if 'S1' in sobol_S1.columns:
                sobol_values = sobol_S1['S1'].values
            else:
                sobol_values = sobol_S1.iloc[:, 0].values

            # RBD-FAST: extract S1 values
            rbd_df = rbd_result.to_df()
            if 'S1' in rbd_df.columns:
                rbd_values = rbd_df['S1'].values
            else:
                rbd_values = rbd_df.iloc[:, 0].values

            # Create comparison dataframe
            comparison_df = pd.DataFrame({
                'Parameter': problem['names'],
                'Sobol': sobol_values,
                'RBD-FAST': rbd_values
            })

            # Create side-by-side bar plot
            fig, ax = plt.subplots(figsize=(12, 6))

            x = np.arange(len(problem['names']))
            width = 0.35

            ax.bar(x - width/2, comparison_df['Sobol'], width, label='Sobol', alpha=0.8, color='skyblue')
            ax.bar(x + width/2, comparison_df['RBD-FAST'], width, label='RBD-FAST', alpha=0.8, color='coral')

            ax.set_xlabel('Parameters', fontsize=12, fontweight='bold')
            ax.set_ylabel('Sensitivity Index (S1)', fontsize=12, fontweight='bold')
            ax.set_title(f'Sobol vs RBD-FAST Comparison - {output_name}', fontsize=14, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(problem['names'], rotation=45, ha='right')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()

            fig.savefig(f"method_comparison/comparison_{output_name}.png", dpi=300, bbox_inches='tight')
            figures[f"{output_name}_comparison"] = fig
            print(f"  ✓ Generated comparison plot for {output_name}")
            plt.close(fig)

            # Save comparison to CSV
            comparison_df.to_csv(f"method_comparison/comparison_{output_name}.csv", index=False)

        except Exception as e:
            print(f"  ✗ Failed to create comparison for {output_name}: {repr(e)}")

    print(f"\n✓ Method comparison complete: {len(figures)} plots generated")
    return figures


# %%------------   INCREMENTAL SAVE FUNCTIONS              ---------------------------------------------------> #
def initialize_results_csv(csv_path, include_index=True):
    """
    Initialize a CSV file for incremental results saving.

    Parameters:
        csv_path: str, path to the CSV file
        include_index: bool, whether to include 'index' column

    Returns:
        None
    """
    columns = ['index', 'w_out', 'LCOH', "LCOH_i"] if include_index else ['w_out', 'LCOH', 'LCOH_i']
    df = pd.DataFrame(columns=columns)
    df.to_csv(csv_path, index=False)
    print(f"Initialized results file: {csv_path}")


def append_results_to_csv(csv_path, results_batch):
    """
    Append a batch of results to CSV file.

    Parameters:
        csv_path: str, path to the CSV file
        results_batch: list of tuples (index, w_out, LCOH)

    Returns:
        None
    """
    df_batch = pd.DataFrame(results_batch, columns=['index', 'w_out', 'LCOH', "LCOH_i"])
    df_batch.to_csv(csv_path, mode='a', header=False, index=False)


def load_completed_indices(csv_path):
    """
    Load indices of already completed evaluations from CSV.

    Parameters:
        csv_path: str, path to the CSV file

    Returns:
        set: set of completed indices
    """
    if not os.path.exists(csv_path):
        return set()

    try:
        df = pd.read_csv(csv_path)
        if 'index' in df.columns and len(df) > 0:
            return set(df['index'].values)
        return set()
    except Exception as e:
        print(f"Warning: Could not load completed indices from {csv_path}: {repr(e)}")
        return set()


def sort_and_finalize_results(csv_path, output_csv_path):
    """
    Sort results by index and save to final output file.

    Parameters:
        csv_path: str, path to the unsorted CSV file
        output_csv_path: str, path to the final sorted CSV file

    Returns:
        None
    """
    df = pd.read_csv(csv_path)
    if 'index' in df.columns:
        df = df.sort_values('index').reset_index(drop=True)
    df.to_csv(output_csv_path, index=False)
    print(f"Finalized and sorted results saved to: {output_csv_path}")


# %%------------   WORKER FUNCTION FOR PARALLEL EXECUTION   ---------------------------------------------------> #
def worker(params_with_idx):
    """
    Worker function to evaluate set_parameters_and_execute for a single parameter combination.
    Uses lazy imports to minimize memory footprint in multiprocessing workers.

    Parameters:
        params_with_idx: tuple (index, params_array)
            - index: int, row index of parameters
            - params_array: array with [k_rock, rho_cp_rock, gradT, pe, mass_flow, depth, l_horiz]

    Returns:
        tuple: (index, w_out, LCOH_ext)
    """
    idx, params = params_with_idx
    try:
        # Ensure project root is in path for imports
        import sys
        import os
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)

        # Lazy imports: only load heavy modules when worker is actually called
        from main_code.well_model.geometry_based_well_models.REELWEEL_model import (
            REELWELLHeatingSection,
            REELWELLGeometry,
            REELWEELBHE,
        )
        from main_code.support.abstract_plant_thermo_point import PlantThermoPoint

        k_rock, rho_cp_rock, gradT, u_inf, mass_flow, depth, l_horiz= params
        w_out, LCOH, lc = set_parameters_and_execute(
            k=k_rock,
            rho_cp=rho_cp_rock,
            grad_t=gradT,
            u_inf=u_inf,
            mass_flow=mass_flow,
            depth=depth,
            l_horiz=l_horiz,
        )
        return idx, w_out, LCOH, lc
    except Exception as e:
        print(f"Worker failed for params index {idx}: {repr(e)}")
        import traceback
        traceback.print_exc()
        # Return default values if computation fails
        return idx, 1e-6, 1e20,1e20
def wrapped_model_parallel(parametres, n_workers=None, csv_path='sobol_results_incremental.csv', save_every=50,
                           timeout_s=600):
    """
    Parallel evaluation of all parameter combinations using multiprocessing with incremental saving.
    WARNING: This may cause memory issues on Windows systems with heavy dependencies.

    Parameters:
        parametres: array of parameter combinations from sobol.sample()
        n_workers: number of parallel workers (default: cpu_count() - 1)
        csv_path: str, path to save incremental results
        save_every: int, save results every N evaluations
        timeout_s: int, max seconds per point; points exceeding it are killed and saved as failed

    Returns:
        tuple: (res_w, res_l) - arrays of results for w_out, and LCOH
    """
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    print(f"Starting parallel evaluation with {n_workers} workers...")
    print(f"Results will be saved every {save_every} iterations to: {csv_path}")

    # Check for existing results to resume from
    completed_indices = load_completed_indices(csv_path)
    if completed_indices:
        print(f"Found {len(completed_indices)} already completed evaluations. Resuming...")
    else:
        # Initialize new results file
        initialize_results_csv(csv_path, include_index=True)

    # Prepare parameter list with indices, skipping completed ones
    params_with_idx = [(idx, p) for idx, p in enumerate(parametres) if idx not in completed_indices]

    if not params_with_idx:
        print("All evaluations already completed! Loading results...")
        df_all = pd.read_csv(csv_path)
        df_all = df_all.sort_values('index').reset_index(drop=True)
        return df_all['w_out'].values, df_all['LCOH'].values

    print(f"Evaluating {len(params_with_idx)} remaining parameter combinations...")

    run_parallel_with_timeout(params_with_idx, n_workers, csv_path,
                              save_every=save_every, timeout_s=timeout_s)

    # Load all results and sort by index
    print(f"\nLoading all results from {csv_path}...")
    df_all = pd.read_csv(csv_path)
    df_all = df_all.sort_values('index').reset_index(drop=True)

    # Extract sorted arrays
    res_w = df_all['w_out'].values

    res_l = df_all['LCOH'].values

    return res_w, res_l


def wrapped_model_sequential(parametres, csv_path='sobol_results_incremental.csv', save_every=50):
    """
    Sequential evaluation of all parameter combinations with incremental saving.

    Parameters:
        parametres: array of parameter combinations from sobol.sample()
        csv_path: str, path to save incremental results
        save_every: int, save results every N evaluations

    Returns:
        tuple: (res_w, res_l) - arrays of results for w_out, and LCOH
    """
    print(f"Starting sequential evaluation of {len(parametres)} parameter combinations...")
    print(f"Results will be saved every {save_every} iterations to: {csv_path}")

    # Check for existing results to resume from
    completed_indices = load_completed_indices(csv_path)
    if completed_indices:
        print(f"Found {len(completed_indices)} already completed evaluations. Resuming...")
    else:
        # Initialize new results file
        initialize_results_csv(csv_path, include_index=True)

    results_buffer = []
    results_dict = {}  # Store results by index for final sorted output

    try:
        for idx, p in enumerate(tqdm(parametres, desc="Evaluating parameter combinations")):
            # Skip if already computed
            if idx in completed_indices:
                continue

            result = worker((idx, p))
            results_buffer.append(result)
            results_dict[result[0]] = result

            # Save batch
            if len(results_buffer) >= save_every:
                append_results_to_csv(csv_path, results_buffer)
                print(f"\n  ✓ Saved {len(results_buffer)} results to {csv_path}")
                results_buffer = []

                # Force garbage collection
                import gc
                gc.collect()

        # Save remaining results
        if results_buffer:
            append_results_to_csv(csv_path, results_buffer)
            print(f"\n  ✓ Saved final {len(results_buffer)} results to {csv_path}")

    except KeyboardInterrupt:
        print("\n\nInterrupted by user!")
        if results_buffer:
            append_results_to_csv(csv_path, results_buffer)
            print(f"  ✓ Saved {len(results_buffer)} pending results before exit")
        raise
    finally:
        # Ensure garbage collection happens
        import gc
        gc.collect()

    # Load all results and sort by index
    print(f"\nLoading all results from {csv_path}...")
    df_all = pd.read_csv(csv_path)
    df_all = df_all.sort_values('index').reset_index(drop=True)

    # Extract sorted arrays
    res_w = df_all['w_out'].values
    res_l = df_all['LCOH'].values

    return res_w, res_l


def wrapped_model_parallel(parametres, n_workers=None, csv_path='sobol_results_incremental.csv', save_every=50,
                           timeout_s=600):
    """
    Parallel evaluation of all parameter combinations using multiprocessing with incremental saving.
    WARNING: This may cause memory issues on Windows systems with heavy dependencies.

    Parameters:
        parametres: array of parameter combinations from sobol.sample()
        n_workers: number of parallel workers (default: cpu_count() - 1)
        csv_path: str, path to save incremental results
        save_every: int, save results every N evaluations
        timeout_s: int, max seconds per point; points exceeding it are killed and saved as failed

    Returns:
        tuple: (res_w, res_l) - arrays of results for w_out, and LCOH
    """
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    print(f"Starting parallel evaluation with {n_workers} workers...")
    print(f"Results will be saved every {save_every} iterations to: {csv_path}")

    # Check for existing results to resume from
    completed_indices = load_completed_indices(csv_path)
    if completed_indices:
        print(f"Found {len(completed_indices)} already completed evaluations. Resuming...")
    else:
        # Initialize new results file
        initialize_results_csv(csv_path, include_index=True)

    # Prepare parameter list with indices, skipping completed ones
    params_with_idx = [(idx, p) for idx, p in enumerate(parametres) if idx not in completed_indices]

    if not params_with_idx:
        print("All evaluations already completed! Loading results...")
        df_all = pd.read_csv(csv_path)
        df_all = df_all.sort_values('index').reset_index(drop=True)
        return df_all['w_out'].values, df_all['LCOH'].values

    print(f"Evaluating {len(params_with_idx)} remaining parameter combinations...")

    run_parallel_with_timeout(params_with_idx, n_workers, csv_path,
                              save_every=save_every, timeout_s=timeout_s)

    # Load all results and sort by index
    print(f"\nLoading all results from {csv_path}...")
    df_all = pd.read_csv(csv_path)
    df_all = df_all.sort_values('index').reset_index(drop=True)

    # Extract sorted arrays
    res_w = df_all['w_out'].values

    res_l = df_all['LCOH'].values

    return res_w, res_l



