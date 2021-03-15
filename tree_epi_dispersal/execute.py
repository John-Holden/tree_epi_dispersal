"""
Run ensembles averaging methods on HPC or local machine.
"""
import datetime
from typing import Union
from tree_epi_dispersal.model_dynamics import run_SIR, run_ADB
from tree_epi_dispersal.ensemble_simulation_helpers import time_print
from tree_epi_dispersal.ensemble_analysis import process_avg_R0_struct
from parameters_and_settings import ParamsAndSetup


def get_avg_R0(rho: float, beta: float) -> list:
    """
    Repeat R0 simulations over a number of `runs' for a given value of rho, beta.
    """
    assert ParamsAndSetup['params'].ensemble_size, 'Error, did not set ensemble-size!'
    ensemble_size = ParamsAndSetup['params'].ensemble_size
    ell = ParamsAndSetup['params'].ell

    ensemble_R0 = []
    for repeat in range(ensemble_size):
        if ParamsAndSetup['setup'].verb >= 2:
            print('Repeat : {}'.format(repeat))
        if ParamsAndSetup['params'].model == 'SIR':
            sim_result = run_SIR(rho, beta, ell)
        elif ParamsAndSetup['params'].model == 'ADB':
            sim_result = run_ADB(rho, beta, ell)
        else:
            raise ValueError(f'Expected model /in [SIR, ADB] found {ParamsAndSetup["params"].model}')
        mean_R0_per_gen = process_avg_R0_struct(sim_result['R0_hist'])

        ensemble_R0.append(mean_R0_per_gen)  # the number of gen-0 secondary infections

    return ensemble_R0


def single_sim(rho: float, beta: float, ell: Union[int, float, tuple], dispersal_type: str = 'gaussian',
               model: str = 'SIR', plot_show: bool = False, plot_freq: int = 10):
    """
    Run a single instance of the dispersal_model."""

    ParamsAndSetup['params'].ell = ell
    ParamsAndSetup['params'].dispersal_model = dispersal_type
    ParamsAndSetup['params'].assert_correct_dispersal()

    if plot_show:
        ParamsAndSetup['setup'].plot = True
        ParamsAndSetup['setup'].show = True
        ParamsAndSetup['setup'].plot_freq = plot_freq

    start = datetime.datetime.now()
    print(f'\n Running {model} | {dispersal_type} @ single_sim')
    print(f'\t beta = {round(beta, 3)}, rho = {round(rho, 3)}')
    if ParamsAndSetup['setup'].verb:
        print(f'\t Model : {ParamsAndSetup["params"].dispersal_model}')

    # run simulation
    if model == 'SIR':
        out = run_SIR(rho, beta, ell)
    elif model == 'ADB':
        ParamsAndSetup['params'].adb_config()
        out = run_ADB(rho, beta, ell)
    else:
        raise ValueError('Wrong dispersal_model input: Implement [SIR, ADB]')

    if 'mortality_ratio' in out:
        print(f'\t Mortality ratio : {out["mortality_ratio"]}')
    if 'R0_hist' in out:
        from tree_epi_dispersal.ensemble_analysis import process_avg_R0_struct
        from tree_epi_dispersal.plot_methods import pltR0
        print(f'R0 hist {out["R0_hist"]}')
        R0_data = process_avg_R0_struct(out['R0_hist'])
        if len(R0_data) > 1:
            pltR0(R0_data)
        else:
            print(f'mean R0 = {R0_data}')

    elapsed = datetime.datetime.now() - start
    print(out['termination'])
    time_print(elapsed.seconds, msg='@ singleSim Done')



