import datetime
import numpy as np
from model import runSim
from timeit import default_timer as timer
from run_model import ModelParamSet, Settings, Metrics

def singleSim(rho, beta):
    """
    Run a single instance of the model.
    """
    from helper_functions import timerPrint, R0_generation_mean

    start = timer()
    print('\n Running @ singleSim...')
    print(f'\t beta = {round(beta, 3)}, rho = {round(rho, 3)}')
    out = runSim(pc=ModelParamSet(rho, beta, alpha=10, L=1000), metrics=Metrics(), settings=Settings())
    elapsed = round(timer() - start, 2)
    print(f'\n@ singleSim DONE | {timerPrint(elapsed)}')
    return out

def R0_analysis(metrics, save=False):
    "plot R0 as a function of generation"
    from helper_functions import  R0_generation_mean
    from plots.plotLib import pltR0
    meanR0_vs_gen = R0_generation_mean(metrics.R0_trace)
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
    from helper_functions import timerPrint
    from ensemble_methods import save_ensemble, save_sim_info, save_sim_out

    if jobId == '1' or jobId == 'local_test':
        save_sim_info(ens_field_names=['R0_trace', 'extinction_time', 'mortality_ratio'],
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
    from ensemble_methods import runR0_ensemble, mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-local-ensemble'
    ens_name = mk_new_dir(ens_name)
    Pspace_iterator(runR0_ensemble, repeats, rhos, betas, ens_name, jobId='local_test')
    return 'Success'

def R0_domain_sensitivity(runs, rho, beta, box_sizes ,Mparams, Mts, sts):
    """
    Run sensitivity analysis on model for different grid sizes.
    """
    import sys
    from collections import defaultdict
    import matplotlib.pyplot as plt
    from helper_functions import timerPrint, R0_generation_mean, avg_multi_dim
    from timeit import default_timer as timer
    R0_gen_ens = defaultdict(list)
    for N in range(runs):
        print('_')
        for iter_, L in enumerate(box_sizes):
            mprms = Mparams(rho=rho, beta=beta, L=L, alpha=5)
            if sts.verbose >= 1:
                print(f'\t Repeat : {N}, Box size : {L}')
                print(f'\t equivalent grid size {L*mprms.alpha/1000}km x {L*mprms.alpha/1000}km')
            if sts.verbose >= 1:
                start = timer()

            [out_mPrm, out_mts] = runSim(pc=mprms, metrics=Mts(), settings=sts)
            meanR0_vs_gen = R0_generation_mean(out_mts.R0_trace)
            if sts.verbose >= 1:
                print('\n\tTime steps elapsed = {}'.format(out_mts.endT))
                print('\tPercolation = {} @ time = {}'.format(out_mts.percolation, out_mts.percT))
                print('\tExtinction = {} @ time = {}'.format(out_mts.extinction, out_mts.extinctionT))
                print('\tMortality Ratio = {} '.format(out_mts.mortality_ratio))
                print('\tTot Number Removed = {} '.format(out_mts.numR[-1]))
                elapsed = round(timer() - start, 2)
                print(f'\t@ singleSim DONE | {timerPrint(elapsed)}\n')
            R0_gen_ens[iter_].append(meanR0_vs_gen)

    for box_size in R0_gen_ens:
        plt.plot(avg_multi_dim(R0_gen_ens[box_size]), label=f' size: {box_sizes[box_size]}, '
                                                            f'grid size:{box_sizes[box_size]*mprms.alpha/1000}km')
    plt.xlabel('Generation', size=20)
    plt.ylabel(r'$\overline{R}_0$', size=25)
    plt.legend()
    plt.show()
