import os
import json
import datetime
from typing import Callable, Union, Iterable

from collections import defaultdict
from tree_epi_dispersal.ensemble_simulation_helpers import save_sim_out, time_print

from parameters_and_settings import PATH_TO_DATA_STORE

def save_output(method:Callable):

    def wrapper(ensemble_name: str, jobId: Union[None, str] = None, **kwargs):



        start = datetime.datetime.now()
        ensemble_output = method(**kwargs)
        end = datetime.datetime.now() - start
        time_print(end)

        if jobId == '1' or jobId == 'local_test':  # write elapsed time to file
            save_sim_out(path_to_ensemble=path_to_ensemble, elapsed_time=timerPrint(elapsed))


        path_to_ensemble = f'{PATH_TO_DATA_STORE}{ensemble_name}/core_output/'
        save_name = f"/core_{jobId}" if jobId else f"/local.json"

        with open(f"{path_to_ensemble}{save_name}",
                  'w') as json_file:  # save ensemble in json struct
            json.dump(ensemble_output, json_file, indent=4)

    return wrapper



@save_output
def parameter_space_iterator(method: Callable, number_samples: int, rhos:Iterable, betas:Iterable) -> "Success":
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations.
    - ensemble_method is a Callable method of the type of ensemble we wish to run e.g. velocity/percolation/time-series
    """

    ensemble_results = defaultdict(dict)
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            all_ens_fields = method(rho, beta, runs=number_samples)
            ensemble_results[f'rho_{rho}_beta_{beta}'] = all_ens_fields

    return ensemble_results