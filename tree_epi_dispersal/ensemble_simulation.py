import json
import datetime
from typing import Callable, Union, Iterable

from collections import defaultdict
from parameters_and_settings import ParamsAndSetup
from tree_epi_dispersal.ensemble_simulation_helpers import save_meta_data, time_print

from parameters_and_settings import PATH_TO_DATA_STORE


def save_output(ensemble_averager: Callable):
    """
    for a given ensemble averaging method, save results and log executed time.
    """
    def wrapper(execute_model: Callable, ensemble_name: str, jobId: Union[None, str] = None):
        """
         Wrapper to save output and display time taken.
        """
        start = datetime.datetime.now()
        ensemble_output = ensemble_averager(execute_model)

        end = datetime.datetime.now() - start
        elapsed = time_print(end.seconds)

        path_to_ensemble = f'{PATH_TO_DATA_STORE}{ensemble_name}'
        save_name = f"core_{jobId}" if jobId else f"/local.json"

        if jobId == '1' or jobId is None:  # write elapsed time to file
            save_meta_data(path_to_ensemble=path_to_ensemble, elapsed_time=elapsed)

        with open(f"{path_to_ensemble}/core_output/{save_name}",
                  'w') as json_file:  # save ensemble in json struct
            json.dump(ensemble_output, json_file, indent=4)

    return wrapper


@save_output
def parameter_space_iterator(execute_model: Callable):
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations.
    - ensemble_method is a Callable method of the type of ensemble we wish to run e.g. velocity/percolation/time-series
    """

    rhos = ParamsAndSetup['params'].rhos
    betas = ParamsAndSetup['params'].betas
    assert len(rhos) > 0 and len(betas) > 0

    ensemble_results = defaultdict(dict)
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            all_ens_fields = execute_model(rho, beta)
            ensemble_results[f'rho_{rho}_beta_{beta}'] = all_ens_fields

    return ensemble_results
