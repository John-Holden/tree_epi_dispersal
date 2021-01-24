import datetime
import numpy as np
from model import runSim
from typing import Type
from timeit import default_timer as timer
from run_simulation import ModelParamSet, Settings, Metrics

def singleSim(rho:float, beta:float, L=1000) -> '([S,I,R], metrics)':
    """
    Run a single instance of the model.
    """
    from helper_methods import timerPrint

    start = timer()
    print('\n Running @ singleSim...')
    print(f'\t beta = {round(beta, 3)}, rho = {round(rho, 3)}')
    out = runSim(pc=ModelParamSet(rho, beta, alpha=5, L=L), metrics=Metrics(), settings=Settings())
    elapsed = round(timer() - start, 2)
    print(f'\n@ singleSim DONE | {timerPrint(elapsed)}')
    return out

def R0_analysis(metrics:Type[Metrics], save=False) -> 'Success':
    "Given metrics class, from singleSim, plot R0 as a function of generation"
    from helper_methods import  R0_generation_mean
    from plots.plotLib import pltR0
    meanR0_vs_gen = R0_generation_mean(metrics.R0_histories)
    print(meanR0_vs_gen, 'mean R0 gen')
    pltR0(meanR0_vs_gen, save)
    print('\n...Time steps elapsed = {}'.format(metrics.endT))
    print('...Percolation = {} @ time = {}'.format(metrics.percolation, metrics.percT))
    print('...Extinction = {} @ time = {}'.format(metrics.extinction, metrics.extinctionT))
    print('...Mortality Ratio = {} '.format(metrics.mortality_ratio))
    print('...Tot Number Removed = {} '.format(metrics.numR[-1]))
    return "Success"


def Pspace_iterator(run_ensemble, N, rhos, betas, ensName=None, jobId=None) -> "Success":
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations
    :param run_ensemble: method, type of ensemble running e.g. velocity/percolation/time-series
    :param rhos: float, av tree densities
    :param betas: float, infectivity constant
    """
    from helper_methods import timerPrint
    from ensemble_averaging_methods import save_ensemble, save_ens_info, save_sim_out

    if jobId == '1' or jobId == 'local_test':
        save_ens_info(ens_field_names=['R0_trace', 'extinction_time', 'mortality_ratio'],
                  rhos=rhos, betas=betas, param_set=ModelParamSet(0, 0), settings=Settings(),
                  ensemble_name=ensName, per_core_repeats=N)

    start = timer()  # Time ensemble averaging process
    R0_results = np.zeros(shape=(len(betas), len(rhos), N))
    extinctionT_results = np.zeros(shape=(len(betas), len(rhos), N))
    mortality_ratio_results = np.zeros(shape=(len(betas), len(rhos), N))
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            all_ens_fields = run_ensemble(rho, beta, runs=N, MPrm=ModelParamSet, Mts=Metrics, Sts=Settings)
            R0_results[i, j] = all_ens_fields[0]
            extinctionT_results[i, j] = all_ens_fields[1]
            mortality_ratio_results[i, j] = all_ens_fields[2]

    save_ensemble(ens_field_names=['R0_trace', 'extinction_time', 'mortality_ratio'],
                 ens_results=[R0_results, extinctionT_results, mortality_ratio_results],
                 ensemble_name=ensName, job_id=jobId)

    elapsed = round(timer() - start, 2)
    if jobId == '1' or jobId == 'local_test':
        save_sim_out(ensemble_name=ensName, elapsed_time=timerPrint(elapsed))
    print('\n@ Ensemble Run DONE | {}'.format(timerPrint(elapsed)))
    return "Success"


def run_lcl_ens(repeats, rhos, betas):
    from ensemble_averaging_methods import runR0_ensemble, mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    Pspace_iterator(runR0_ensemble, repeats, rhos, betas, ens_name, jobId='local_test')
    return 'Success'


def run_lcl_R0_sensitivity(repeats, rho, beta, alpha, box_sizes):
    from ensemble_averaging_methods import mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    R0_domain_sensitivity(runs=repeats, rho=rho, beta=beta, alpha=alpha, box_sizes=box_sizes, jobId='local-run', ens_name=ens_name)



def R0_domain_sensitivity(runs:int, rho:float, beta:float, alpha:int, box_sizes:list, jobId:str, ens_name=None):
    """
    Run sensitivity analysis on model for different grid sizes.
    """
    import os
    import json
    from collections import defaultdict
    from ensemble_averaging_methods import save_sim_out, save_ens_info
    from helper_methods import timerPrint, R0_generation_mean
    from timeit import default_timer as timer
    R0_gen_ens = defaultdict(list)
    start = timer()
    if jobId == '1' or jobId == 'local-run':  # if hpc simulation or local test
        save_ens_info(ens_field_names=['R0_histories'], box_sizes=box_sizes,
                      rhos=rho, betas=beta, param_set=ModelParamSet(rho=rho, beta=beta, alpha=alpha), settings=Settings(),
                      ensemble_name=ens_name, per_core_repeats=runs)
    sim_settings = Settings()
    for N in range(runs):
        for iter_, L in enumerate(box_sizes):
            model_params = ModelParamSet(rho=rho, beta=beta, L=L, alpha=alpha)
            if sim_settings.verbose >= 1:
                print(f'\t Repeat : {N}, Box size : {L}')
                print(f'\t equivalent grid size {L*model_params.alpha/1000}km x {L*model_params.alpha/1000}km')
            out_mts = runSim(pc=model_params, metrics=Metrics(), settings=sim_settings)[1]
            meanR0_vs_gen = R0_generation_mean(out_mts.R0_histories)
            R0_gen_ens[L].append(meanR0_vs_gen)

    elapsed = round(timer() - start, 2)
    print(f'\t@ singleSim DONE | {timerPrint(elapsed)}\n')

    if jobId == '1' or jobId == 'local-run':
        save_sim_out(ensemble_name=ens_name, elapsed_time=timerPrint(elapsed))

    with open(f"{os.getcwd()}/ensemble_dat/{ens_name}/R0_histories/core_{jobId}.json", 'w') as fp:
        json.dump(R0_gen_ens, fp, indent=4)  # save as json struct

    return 'Success'

