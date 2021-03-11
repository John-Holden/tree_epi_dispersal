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
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-hpc-R0-vs-rho'
    mk_new_dir(ens_name, job_id)

    ParamsAndSetup['setup'].plot = False
    ParamsAndSetup['setup'].verb = 0
    ParamsAndSetup['setup'].max_generation_bcd = 3

    Ensemble = ParamsAndSetup['ensemble']
    ens_conf = Ensemble('part')
    ParamsAndSetup['params'].rhos = ens_conf.rhos
    ParamsAndSetup['params'].betas = ens_conf.betas
    ParamsAndSetup['params'].L = 1000
    ParamsAndSetup['params'].ell = 195
    ParamsAndSetup['params'].model = 'gaussian'
    ParamsAndSetup['params'].ensemble_size = 1
    ParamsAndSetup['params'].update_epi_c()
    ParamsAndSetup['params'].assert_config()

    beta_rho_iterator(execute_model=get_avg_R0, ensemble_name=ens_name, job_id=job_id)


def local_mode():
    """
    Run ensemble simulation on local machine
    """
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    mk_new_dir(ens_name)

    ParamsAndSetup['setup'].max_generation_bcd = 3
    ParamsAndSetup['setup'].plot = False
    ParamsAndSetup['setup'].verb = 1

    ParamsAndSetup['params'].rhos = [0.01, 0.005]
    ParamsAndSetup['params'].betas = [0.05]
    ParamsAndSetup['params'].ensemble_mode = True
    ParamsAndSetup['params'].ensemble_size = 2

    beta_rho_iterator(execute_model=get_avg_R0, ensemble_name=ens_name)


if __name__ == '__main__':
    hpc_mode_ = len(sys.argv) > 1
    if 'HPC_MODE' in os.environ:
        assert len(sys.argv) > 1
        hpc_mode()
    else:
        local_mode()



