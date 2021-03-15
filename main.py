import os
import sys
import datetime

from tree_epi_dispersal.ensemble_simulation_helpers import mk_new_dir
from tree_epi_dispersal.execute import get_avg_R0
from tree_epi_dispersal.ensemble_simulation import beta_rho_iterator
from parameters_and_settings import ParamsAndSetup


def hpc_mode():
    """
    Run simulation via hpc task array
    """

    job_id = sys.argv[1:][0]
    ParamsAndSetup['setup'].ensemble_config()

    Ensemble = ParamsAndSetup['ensemble']
    ens_conf = Ensemble('part')
    ParamsAndSetup['params'].rhos = ens_conf.rhos
    ParamsAndSetup['params'].betas = [0.000050, 0.00015]

    ParamsAndSetup['params'].ell = 195
    ParamsAndSetup['params'].ensemble_size = 5
    ParamsAndSetup['params'].dispersal_model = 'gaussian'
    ParamsAndSetup['params'].ensemble_mode = True
    ParamsAndSetup['params'].adb_config()  # set config to model ash dieback
    ParamsAndSetup['params'].assert_config()

    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = f'{date}-hpc-R0-vs-rho-{ParamsAndSetup["params"].model}'
    mk_new_dir(ens_name, job_id)

    beta_rho_iterator(execute_model=get_avg_R0, ensemble_name=ens_name, job_id=job_id)


def local_mode():
    """
    Run ensemble simulation on local machine
    """

    Ensemble = ParamsAndSetup['ensemble']
    ens_conf = Ensemble('part')

    ParamsAndSetup['params'].rhos = ens_conf.rhos
    ParamsAndSetup['params'].betas = [0.0001, 0.001]
    ParamsAndSetup['params'].ensemble_mode = True
    ParamsAndSetup['params'].ensemble_size = 10

    ParamsAndSetup['params'].adb_config()  # set config to model ash dieback

    ParamsAndSetup['setup'].ensemble_config()  # set config to ensemble mode
    ParamsAndSetup['setup'].verb = 1

    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = f'{date}-local-ensemble-{ParamsAndSetup["params"].model}'
    mk_new_dir(ens_name, job_id='local')

    beta_rho_iterator(execute_model=get_avg_R0, ensemble_name=ens_name)


if __name__ == '__main__':
    hpc_mode_ = len(sys.argv) > 1
    if 'HPC_MODE' in os.environ:
        assert len(sys.argv) > 1
        hpc_mode()
    else:
        local_mode()



