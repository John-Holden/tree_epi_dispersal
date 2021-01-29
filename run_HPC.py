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
ensemble_name = date+'-hpc-R0-vs-rho-test-two-cluster-sizes'
if job_id == '1':
    ensemble_name = mk_new_dir(ensemble_name)
# - make directories & run phase
# ----------------------------------------------------- #
assert Settings().plot == False

rho_small = np.arange(0.0, 0.031, 0.001)
rho_large = np.arange(0.04, 0.101, 0.01)
rhos = np.hstack([rho_small, rho_large])
betas = np.array([4e-05, 6.e-05]) # are cluster sizes the same for two values of beta ?

result = parameter_space_iterator(ensemble_method=run_R0_ensemble, ensName=ensemble_name, jobId=job_id,
                                  N=10, rhos=rhos, betas=betas)


