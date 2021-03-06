"""
Run ensembles averaging methods on HPC or local machine.
"""
from typing import Union
from tree_epi_dispersal.model_dynamics import run_simulation
from timeit import default_timer as timer
from parameters_and_settings import ParamsAndSetup


def get_avg_R0(rho: float, beta: float) -> dict:
    """
    Repeat R0 simulations over a number of `runs' for a given value of rho, beta.
    """
    from tree_epi_dispersal.model_dynamics_helpers import R0_generation_mean
    ensemble_R0 = []
    ensemble_size = ParamsAndSetup['params'].ensemble_size
    ell = ParamsAndSetup['params'].ell
    for repeat in range(ensemble_size):
        if ParamsAndSetup['setup'].verb >= 1:
            print('Repeat : {}'.format(repeat))

        sim_result = run_simulation(rho, beta, ell)
        assert 0
        R0_histories = R0_generation_mean(out_metrics._R0_histories)
        ensemble_R0.append(R0_histories)  # the number of gen-0 secondary infections
    return {'mean_R0_vs_gen_core_ensemble': ensemble_R0}


def R0_analysis(save=False) -> 'Success':
    "Given metrics class, from singleSim, plot R0 as a function of generation"
    from tree_epi_dispersal.model_dynamics_helpers import  R0_generation_mean
    from tree_epi_dispersal.plot_methods import pltR0
    Metrics = ParamsAndSetup['metrics']
    metrics = Metrics()
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
    run_simulation(pc=ModelParamSet(rho, beta, alpha=5, L=L), metrics=Metrics())
    elapsed = round(timer() - start, 2)
    print(f'\n@ singleSim DONE | {timerPrint(elapsed)}')


