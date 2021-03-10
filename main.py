import os
import sys
import datetime
import numpy as np
from parameters_and_settings import Settings

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
    ensemble_name = date + '-hpc-R0-vs-rho'
    if job_id == '1':
        mk_new_dir(ensemble_name)

    assert Settings().plot == False

    rho_small = np.arange(0.0, 0.031, 0.0025)
    rho_large = np.arange(0.04, 0.101, 0.01)
    rhos = np.hstack([rho_small, rho_large])  # len 20
    betas = np.arange(0.00000, 0.00032, 0.00002)  # len 16

    # parameter_space_iterator(method=g, ensemble_name=ensemble_name, jobId=job_id,
    #                          number_samples=2, rhos=rhos, betas=betas)


def local_mode():
    """
    Run ensemble simulation on local machine
    """
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    mk_new_dir(ens_name)

    ParamsAndSetup['params'].rhos = [0.99, 0.005]
    ParamsAndSetup['params'].betas = [0.05]
    ParamsAndSetup['params'].ensemble_mode = True
    ParamsAndSetup['params'].ensemble_size = 2

    beta_rho_iterator(execute_model=get_avg_R0, ensemble_name=ens_name)


if __name__ == '__main__':
    hpc_mode_ = len(sys.argv) > 1
    if 'HPC_MODE' in os.environ:
        assert len(sys.argv)>1
        hpc_mode()
    else:
        local_mode()



