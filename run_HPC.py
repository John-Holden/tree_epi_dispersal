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
result = parameter_space_iterator(ensemble_method=run_R0_ensemble, ensName=ensemble_name, jobId=job_id,
                                  N=2, rhos=np.arange(0.002, 0.102, 0.002), betas=np.arange(0.0001, 0.0011, 0.0001))

