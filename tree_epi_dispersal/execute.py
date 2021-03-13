"""
Run ensembles averaging methods on HPC or local machine.
"""
import datetime
from typing import Union
from tree_epi_dispersal.model_dynamics import run_SIR
from tree_epi_dispersal.ensemble_simulation_helpers import time_print
from parameters_and_settings import ParamsAndSetup


def get_avg_R0(rho: float, beta: float) -> list:
    """
    Repeat R0 simulations over a number of `runs' for a given value of rho, beta.
    """
    from tree_epi_dispersal.model_dynamics_helpers import R0_generation_mean
    ensemble_R0 = []
    assert ParamsAndSetup['params'].ensemble_size, 'Error, did not set ensemble-size!'
    ensemble_size = ParamsAndSetup['params'].ensemble_size
    ell = ParamsAndSetup['params'].ell
    for repeat in range(ensemble_size):
        if ParamsAndSetup['setup'].verb == 2:
            print('Repeat : {}'.format(repeat))
        sim_result = run_SIR(rho, beta, ell)
        mean_R0_per_gen = R0_generation_mean(sim_result['R0_hist'])
        ensemble_R0.append(mean_R0_per_gen)  # the number of gen-0 secondary infections

    return ensemble_R0


def single_sim(rho: float, beta: float, ell: Union[int, float, tuple], model: str = 'gaussian',
               plot_show: bool = False, plot_freq: int = 10):
    """
    Run a single instance of the model."""

    ParamsAndSetup['params'].ell = ell
    ParamsAndSetup['params'].model = model
    ParamsAndSetup['params'].assert_correct_dispersal()

    if plot_show:
        ParamsAndSetup['setup'].plot = True
        ParamsAndSetup['setup'].show = True
        ParamsAndSetup['setup'].plot_freq = plot_freq

    start = datetime.datetime.now()
    print('\n Running @ singleSim...')
    print(f'\t beta = {round(beta, 3)}, rho = {round(rho, 3)}')
    if ParamsAndSetup['setup'].verb:
        print(f'\t Model : {ParamsAndSetup["params"].model}')
    out = run_SIR(rho, beta, ell)

    if 'mortality_ratio' in out:
        print(f'\t Mortality ratio : {out["mortality_ratio"]}')
    if 'R0_hist' in out:
        from tree_epi_dispersal.ensemble_analysis import process_avg_R0_struct
        from tree_epi_dispersal.plot_methods import pltR0

        R0_data = process_avg_R0_struct(out['R0_hist'])
        pltR0(R0_data)

    elapsed = datetime.datetime.now() - start
    print(out['termination'])
    time_print(elapsed.seconds, msg='@ singleSim Done')



