"""
Run ensembles averaging methods on HPC or local machine.
"""
from tree_epi_dispersal.model_dynamics import runSim
from typing import Type
from timeit import default_timer as timer
from parameters_and_settings import ModelParamSet, Settings, Metrics


def run_R0_ensemble(rho:float, beta:float, runs:int) -> dict:
    """
    Repeat R0 simulations over a number of `runs' for a given value of rho, beta.
    """
    from tree_epi_dispersal.model_dynamics_helpers import R0_generation_mean
    ensemble_R0 = []
    for N in range(runs):
        if Settings.verbose:
            print('Repeat : {}'.format(N))
        out_metrics = runSim(pc=ModelParamSet(rho, beta), metrics=Metrics())[1]
        R0_histories = R0_generation_mean(out_metrics.R0_histories)
        ensemble_R0.append(R0_histories)  # the number of gen-0 secondary infections
    return {'mean_R0_vs_gen_core_ensemble': ensemble_R0}


def R0_analysis(metrics:Type[Metrics], save=False) -> 'Success':
    "Given metrics class, from singleSim, plot R0 as a function of generation"
    from tree_epi_dispersal.model_dynamics_helpers import  R0_generation_mean
    from tree_epi_dispersal.plot_methods import pltR0

    meanR0_vs_gen = R0_generation_mean(metrics.R0_histories)
    print(meanR0_vs_gen, 'mean R0 gen')
    pltR0(meanR0_vs_gen, save)
    print('\n...Time steps elapsed = {}'.format(metrics.endT))
    print('...Percolation = {} @ time = {}'.format(metrics.percolation, metrics.percT))
    print('...Extinction = {} @ time = {}'.format(metrics.extinction, metrics.extinctionT))
    print('...Mortality Ratio = {} '.format(metrics.mortality_ratio))
    print('...Tot Number Removed = {} '.format(metrics.numR[-1]))
    return "Success"


def singleSim(rho:float, beta:float, L=1000) -> '([S,I,R], metrics)':
    """
    Run a single instance of the model.
    """
    from tree_epi_dispersal.model_dynamics_helpers import timerPrint

    start = timer()
    print('\n Running @ singleSim...')
    print(f'\t beta = {round(beta, 3)}, rho = {round(rho, 3)}')
    runSim(pc=ModelParamSet(rho, beta, alpha=5, L=L), metrics=Metrics())
    elapsed = round(timer() - start, 2)
    print(f'\n@ singleSim DONE | {timerPrint(elapsed)}')


