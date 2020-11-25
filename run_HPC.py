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
ens_name = date+'-hpc-vel-ens'
if job_id == '1':
    ens_name = mk_new_dir(ens_name)
# - make directories & run phase
# ----------------------------------------------------- #
rhos = np.linspace(0.0, 0.03, 3)
betas = np.linspace(0.003, 0.006, 2)
repeats = 2
result = Pspace_iterator(runR0_ensemble, repeats, rhos, betas, ens_name, job_id)
print('...'+result)
sys.exit()
