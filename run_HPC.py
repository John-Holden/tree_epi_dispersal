import sys
import datetime
from runner_methods import R0_domain_sensitivity
from ensemble_averaging_methods import mk_new_dir

# - Input variables & setup/save directories
# ----------------------------------------------------- #
job_id = sys.argv[1:][0]
date = datetime.datetime.today().strftime('%Y-%m-%d')
ens_name = date+'-hpc-R0-generation'
if job_id == '1':
    ens_name = mk_new_dir(ens_name)
# - make directories & run phase
# ----------------------------------------------------- #
result = R0_domain_sensitivity(runs=5, rho=0.01, beta=.0001, box_sizes=[250, 500, 750, 1000],
                               jobId=job_id,ens_name=ens_name)
print(result)
