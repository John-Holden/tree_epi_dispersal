import json
import datetime
from typing import Callable, Union

from collections import defaultdict
from parameters_and_settings import ParamsAndSetup, PATH_TO_DATA_STORE
from tree_epi_dispersal.ensemble_simulation_helpers import save_meta_data, time_print, write_time


def save_output(ensemble_averager: Callable):
    """
    for a given ensemble averaging method, save results and log executed time.
    """
   
    def wrapper(execute_model: Callable, ensemble_name: str, job_id: Union[None, str] = None):
        """
         Wrapper to save output and display time taken.
        """
        start = datetime.datetime.now()
        path_to_ensemble = f'{PATH_TO_DATA_STORE}{ensemble_name}'
        save_name = f"core_{job_id}.json" if job_id else f"/local.json"

        save_meta_data(path_to_ensemble, job_id)
        ensemble_output = ensemble_averager(execute_model)
        end = datetime.datetime.now() - start
        elapsed = time_print(end.seconds)

        write_time(path_to_ensemble, elapsed, job_id)

        with open(f"{path_to_ensemble}/core_output/{save_name}",
                  'w') as json_file:  # save ensemble in json struct
            json.dump(ensemble_output, json_file, indent=4)

    return wrapper


@save_output
def beta_rho_iterator(execute_model: Callable):
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations.
    - ensemble_method is a Callable method of the type of ensemble we wish to run e.g. velocity/percolation/time-series
    """
    rhos = ParamsAndSetup['params'].rhos
    betas = ParamsAndSetup['params'].betas
    assert len(rhos) > 0 and len(betas) > 0
    c = 0
    freq = (len(rhos) * len(betas) )/ 10
    ensemble_results = defaultdict(dict)
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            if ParamsAndSetup['setup'].verb == 1 and c % freq == 0:
                print(f'\t i : {i} / {len(betas)}, j : {j} / {len(rhos)}')

            all_ens_fields = execute_model(rho, beta)
            ensemble_results[f'rho_{rho}_beta_{beta}'] = all_ens_fields

    return ensemble_results
