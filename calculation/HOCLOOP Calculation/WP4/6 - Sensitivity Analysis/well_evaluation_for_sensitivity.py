"""
Helper script to complete any remaining Sobol evaluations after a crash.
This script will load your existing results and complete only the missing evaluations.
"""

print("Starting well_evaluation_for_sensitivity.py...")
print("Importing libraries...")

import pandas as pd
from tqdm import tqdm
tqdm.monitor_interval = 0
import numpy as np
# Import your model function and CSV helpers
import sys
import os
from scipy.stats import qmc

# Add project root to Python path so main_code can be imported
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"Added project root to path: {project_root}")

sys.path.append(os.path.dirname(__file__))

# Import necessary functions from GSA_Library
from main_code.support.GSA_Library import (
    worker,
    append_results_to_csv,
    load_completed_indices,
    initialize_results_csv,
    run_parallel_with_timeout
)


def complete_missing_evaluations(
    csv_params_file='sobol_parameters_20260225_0106.csv',
    csv_results_file='sobol_results_incremental_20260225_0106.csv',
    save_every=2,
    invert_array=False,
    use_parallel=True,
    n_workers=None,
    timeout_s=600
):
    """
    Complete missing evaluations from a crashed Sobol run using parallel processing.

    Parameters:
        csv_params_file: CSV file with parameter combinations
        csv_results_file: CSV file with partial results
        save_every: Save results every N evaluations
        invert_array: bool, whether to invert the order of missing indices
        use_parallel: bool, whether to use parallel processing (default: True)
        n_workers: int, number of parallel workers (default: cpu_count() - 1)
        timeout_s: max seconds per point (parallel mode); points exceeding it are killed
            and saved as failed (1e-6, 1e20, 1e20), their index logged to <results>_timeouts.txt
    """

    from multiprocessing import cpu_count
    import gc

    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    print("COMPLETING MISSING EVALUATIONS")
    print(f"Parallel mode: {use_parallel}, Workers: {n_workers if use_parallel else 'N/A'}")

    # Load parameters
    print(f"\nLoading parameters from: {csv_params_file}")
    df_params = pd.read_csv(csv_params_file)
    parametres = df_params.values
    print(f"  ✓ Loaded {len(parametres)} total parameter combinations")

    # Load completed indices
    print(f"\nLoading completed results from: {csv_results_file}")
    completed_indices = load_completed_indices(csv_results_file)
    print(f"  ✓ Found {len(completed_indices)} completed evaluations")

    # Initialize CSV if starting fresh
    if not completed_indices and not os.path.exists(csv_results_file):
        initialize_results_csv(csv_results_file, include_index=True)

    # Find missing indices
    all_indices = set(range(len(parametres)))
    missing_indices = sorted(list(all_indices - completed_indices))
    if invert_array:
        missing_indices = missing_indices[::-1]

    if not missing_indices:
        print("\n✓ All evaluations are complete! Nothing to do.")
        return

    print(f"\n⚠ Found {len(missing_indices)} missing evaluations:")
    print(f"  Missing indices: {missing_indices[:20]}{'...' if len(missing_indices) > 20 else ''}")

    # Prepare parameter list with indices for processing
    params_with_idx = [(idx, parametres[idx]) for idx in missing_indices]

    if use_parallel:
        # ... PARALLEL EVALUATION ...
        print(f"\nCompleting missing evaluations (PARALLEL with {n_workers} workers)...")

        print(f"  Timeout per point: {timeout_s}s")
        run_parallel_with_timeout(
            params_with_idx, n_workers, csv_results_file,
            save_every=save_every, timeout_s=timeout_s
        )

    else:
        # ... SEQUENTIAL EVALUATION ...
        print(f"\nCompleting missing evaluations (SEQUENTIAL)...")
        results_buffer = []

        for i, (idx, p) in enumerate(tqdm(params_with_idx, desc="Evaluating parameter combinations"), 1):
            try:
                result = worker((idx, p))
                results_buffer.append(result)

                # Save batch
                if len(results_buffer) >= save_every:
                    append_results_to_csv(csv_results_file, results_buffer)
                    print(f"  ✓ Saved {len(results_buffer)} results")
                    results_buffer = []
                    gc.collect()

            except Exception as e:
                print(f"  ✗ Failed to evaluate index {idx}: {repr(e)}")
                # Keep index alignment by saving defaults on failure
                results_buffer.append((idx, 1e-6, 1e20, 1e20))

        # Save remaining results
        if results_buffer:
            append_results_to_csv(csv_results_file, results_buffer)
            print(f"\n  ✓ Saved final {len(results_buffer)} results")

        gc.collect()

    # Verify completion
    final_completed = load_completed_indices(csv_results_file)
    print(f"\n✓ Completed! Total evaluations: {len(final_completed)}/{len(parametres)}")

    if len(final_completed) < len(parametres):
        still_missing = sorted(list(all_indices - final_completed))
        print(f"⚠ Still missing {len(still_missing)} evaluations: {still_missing[:20]}")
    else:
        print("✓ ALL EVALUATIONS COMPLETE!")


if __name__ == '__main__':

    # Preparing Sampled Parameters using LHS

    CSV_PARAMS = 'parameters_lhc.csv'
    N_LHS = 10000
    LHS_SEED = 42

    COLUMNS = ["Rock Properties", "Rock Properties", 'Target Temperature', "Rock Properties", "Mass Flow",
                   'Target Temperature', ' Horizontal Lengh']
    LO = np.array([1.576534964, 2.08e6, 0.021, -3, 2.5, 2000., 0.])
    HI = np.array([5.259210182, 4196806.26, 0.12, 1, 12.5, 5000., 5000.])

    sampler = qmc.LatinHypercube(d=7, optimization="random-cd", seed=LHS_SEED)
    unit = sampler.random(N_LHS)
    X = qmc.scale(unit, LO, HI)
    df = pd.DataFrame(X, columns=COLUMNS)
    df.to_csv(CSV_PARAMS, index=False)


    CSV_RESULTS = 'results_incremental_lhc.csv'

    # Configuration
    USE_PARALLEL = True  # Set to False for sequential processing
    NUM_WORKERS = None   # None = cpu_count() - 1; set to specific number to limit (e.g., 4, 8)
    SAVE_EVERY = 1      # Save results every N evaluations
    TIMEOUT_S = 600     # Max seconds per point (parallel mode); tune on typical convergence times

    print(f"\nWill complete missing evaluations for:")
    print(f"  Parameters: {CSV_PARAMS}")
    print(f"  Results: {CSV_RESULTS}")

    print(f"  Mode: {'PARALLEL' if USE_PARALLEL else 'SEQUENTIAL'}")
    if USE_PARALLEL and NUM_WORKERS:
        print(f"  Workers: {NUM_WORKERS}")
    print(f"  Save every: {SAVE_EVERY} evaluations")
    if USE_PARALLEL:
        print(f"  Timeout per point: {TIMEOUT_S}s")

    input("\nPress ENTER to continue or Ctrl+C to cancel...")

    complete_missing_evaluations(
        CSV_PARAMS,
        CSV_RESULTS,
        save_every=SAVE_EVERY,
        invert_array=True,
        use_parallel=USE_PARALLEL,
        n_workers=NUM_WORKERS,
        timeout_s=TIMEOUT_S
    )
