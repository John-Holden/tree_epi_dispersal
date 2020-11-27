import sys
import datetime
import os, sys
import numpy as np
from run_model import Pspace_iterator
from ensemble_methods import runR0_ensemble, mk_new_dir

# - Input variables & setup/save directories
# ----------------------------------------------------- #
job_id = sys.argv[1:][0]
date = datetime.datetime.today().strftime('%Y-%m-%d')
ens_name = date+'-hpc-ensemble'
if job_id == '1':
    ens_name = mk_new_dir(ens_name)
# - make directories & run phase
# ----------------------------------------------------- #
rhos = np.linspace(0.00, 0.10, 51)
betas = np.linspace(0.000015, 0.000020, 2)
repeats = 1
result = Pspace_iterator(runR0_ensemble, repeats, rhos, betas, ens_name, job_id)
print('...'+result)
sys.exit()
