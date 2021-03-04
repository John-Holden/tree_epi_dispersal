import os
import sys
import datetime
import numpy as np
from PARAMETERS_AND_SETUP import Settings

from ensemble_averaging_methods import mk_new_dir, run_R0_ensemble
from runner_methods import parameter_space_iterator, R0_domain_sensitivity


def hpc_mode():
    """
    Run simulation on hpc task array
    """
    print(os.environ)
    print(os.environ['HPC_MODE'])

    assert 0
    job_id = sys.argv[1:][0]
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ensemble_name = date + '-hpc-R0-vs-rho'
    if job_id == '1':
        mk_new_dir(ensemble_name)

    assert Settings().plot == False

    rho_small = np.arange(0.0, 0.031, 0.0025)
    rho_large = np.arange(0.04, 0.101, 0.01)
    rhos = np.hstack([rho_small, rho_large])  #  len 20
    betas = np.arange(0.00000, 0.00032, 0.00002)  # len 16

    parameter_space_iterator(ensemble_method=run_R0_ensemble, ensName=ensemble_name, jobId=job_id,
                             N=2, rhos=rhos, betas=betas)


def local_mode():
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    R0_domain_sensitivity(runs=1, rho=0.01, beta=0.001, alpha=5, box_sizes=[1000, 1500], jobId='local-run',
                          ens_name=ens_name)


if __name__ == '__main__':

    hpc_mode = len(sys.argv) > 1
    if hpc_mode:
        hpc_mode()

    else:
        local_mode()



