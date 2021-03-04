import json
import datetime
import numpy as np
from typing import Callable, Union

def parameter_space_iterator(ensemble_method: Callable, N: int,
                             rhos:Union[np.ndarray, list, tuple], betas:Union[np.ndarray, list, tuple],
                             ensName: str, jobId: str) -> "Success":
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations.
    - ensemble_method is a Callable method of the type of ensemble we wish to run e.g. velocity/percolation/time-series
    """
    from collections import defaultdict
    from tree_epi_dispersal.model_dynamics_helpers import timerPrint
    from tree_epi_dispersal.ensemble_simulation_helpers import save_ens_info, save_sim_out

    path_to_ensemble = f'{os.getcwd()}/ensemble_dat/{ensName}'
    if jobId == '1' or jobId == 'local_test':  # write parameters to file before ensemble -- assurance
        save_ens_info(ens_field_names=['R0_Histories'], rhos=rhos, betas=betas, path_to_ensemble=path_to_ensemble,
                      per_core_repeats=N, box_sizes=None)

    start = timer()  # Time ensemble averaging process
    ensemble_results = defaultdict(dict)
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            all_ens_fields = ensemble_method(rho, beta, runs=N)
            ensemble_results[f'rho_{rho}_beta_{beta}'] = all_ens_fields

    elapsed = round(timer() - start, 2)
    if jobId == '1' or jobId == 'local_test':  # write elapsed time to file
        save_sim_out(path_to_ensemble=path_to_ensemble, elapsed_time=timerPrint(elapsed))

    with open(f"{path_to_ensemble}/core_output/core_{jobId}.json", 'w') as json_file:  # save ensemble in json struct
        json.dump(ensemble_results, json_file, indent=4)

    print('\n@ Ensemble Run DONE | {}'.format(timerPrint(elapsed)))
    return "Success"