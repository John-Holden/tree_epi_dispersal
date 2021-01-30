import sys
import datetime
import numpy as np
from PARAMETERS_AND_SETUP import Settings
from runner_methods import parameter_space_iterator
from ensemble_averaging_methods import run_R0_ensemble
from ensemble_averaging_methods import mk_new_dir

# - Input variables & setup/save directories
# ----------------------------------------------------- #
job_id = sys.argv[1:][0]
date = datetime.datetime.today().strftime('%Y-%m-%d')
ensemble_name = date+'-hpc-R0-vs-rho'
if job_id == '1':
    ensemble_name = mk_new_dir(ensemble_name)
# - make directories & run phase
# ----------------------------------------------------- #
assert Settings().plot == False

rho_small = np.arange(0.0, 0.031, 0.0025)
rho_large = np.arange(0.04, 0.101, 0.01)
rhos = np.hstack([rho_small, rho_large])  # len 20
betas = np.arange(0.00000, 0.00032, 0.00002)  # len 16

result = parameter_space_iterator(ensemble_method=run_R0_ensemble, ensName=ensemble_name, jobId=job_id,
                                  N=2, rhos=rhos, betas=betas)


