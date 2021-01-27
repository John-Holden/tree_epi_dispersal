"""
Script to run local simulations.
"""
from typing import Union
import numpy as np
from runner_methods import parameter_space_iterator, R0_domain_sensitivity
import datetime


def run_lcl_ens(repeats: int, rhos: Union[np.ndarray, list, tuple], betas:Union[np.ndarray, list, tuple]) -> 'Success':
    """
    Run local ensemble average over beta and rho.
    """
    from ensemble_averaging_methods import run_R0_ensemble, mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    return parameter_space_iterator(run_R0_ensemble, repeats, rhos, betas, ens_name, jobId='local_test')


def run_lcl_R0_sensitivity(repeats: int, rho:float, beta:float, alpha:Union[float, int],
                           box_sizes:Union[np.ndarray, list, tuple]):
    """
    Run local domain sensitivity over different values of box size L.
    """
    from ensemble_averaging_methods import mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    return R0_domain_sensitivity(runs=repeats, rho=rho, beta=beta, alpha=alpha, box_sizes=box_sizes, jobId='local-run',
                          ens_name=ens_name)

if __name__ == '__main__':
    from runner_methods import singleSim, R0_analysis
    R0_analysis(singleSim(rho=0.05, beta=0.0003)[1])
