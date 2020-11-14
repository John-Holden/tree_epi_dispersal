import sys
import datetime
import os, sys
import numpy as np
from run_model import getPspace
from ensemble_methods import runVel_ensemble, mk_new_dir

# - Input variables & setup/save directories
# ----------------------------------------------------- #
job_id = sys.argv[1:][0]
date = datetime.datetime.today().strftime('%Y-%m-%d')
ens_name = date+'-hpc-vel-ens'
ens_name = mk_new_dir(ens_name)
# - make directories & run phase
# ----------------------------------------------------- #
rhos = np.linspace(0.01, 0.05, 5)
betas = [0.020]
repeats = 1
result = getPspace(runVel_ensemble, repeats, rhos, betas, ens_name, job_id)
print('...'+result)
sys.exit()
