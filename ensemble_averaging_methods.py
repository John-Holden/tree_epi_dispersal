import model_dynamics
from typing import Union, Type
import numpy as np
perc = lambda p: 1 if p else 0
from PARAMETERS_AND_SETUP import Metrics, ModelParamSet, Settings

def mk_new_dir(name: str) -> str:
    "Save new directory to file."
    import os
    import sys
    if os.path.exists(f'{os.getcwd()}/ensemble_dat/{name}'):
        sys.exit(f'ERROR DUPLICATE DIRECTORY{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/core_output')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/info/')
    return name


def save_ens_info(ens_field_names: list, rhos: Union[np.ndarray, float], betas:Union[np.ndarray, float],
                  param_set, settings, path_to_ensemble:str, per_core_repeats:int, box_sizes=None):
    """
    Save simulation ensemble to file at the beginning of the simulation.
    """
    from datetime import datetime

    save_info = {"\t ensemble averaged metrics: ": ens_field_names,
                 "\t start time: ": datetime.now(),
                 "\t model: ": param_set.model,
                 "\t alpha: ": str(param_set.alpha) + '(m)',
                 "\t ell: ": str(param_set.ell) + '(m)',
                 "\t core repeats: ": per_core_repeats,
                 "\t initial epicenter radius: ": str(param_set.r),
                 "\t percolation boundary: ": settings.boundary}

    try:  # save iterable parameter types.
        iter(rhos)
        np.save(f"{path_to_ensemble}/info/rhos", rhos)
        save_info['\t rhos'] = f'{rhos[0]}, {rhos[-1]} | {len(rhos)}'
    except TypeError:
        save_info["\t rho"] = f'{round(rhos, 5)}'
    try:
        iter(betas)
        np.save(f"{path_to_ensemble}/info/betas", betas)
        save_info["\t betas"] = f'{betas[0]}, {betas[-1]} | {len(betas)}'
    except TypeError:
        save_info["\t beta"] = f'{round(betas, 7)}'
    try:
        iter(box_sizes)
        np.save(f"{path_to_ensemble}/info/box_sizes", box_sizes)
        save_info["\t domain size "] = [L for L in box_sizes]
        save_info["\t modelled size"] = [f'{L * param_set.alpha / 1000} x {L * param_set.alpha / 1000} (km)' for L in
                                         box_sizes]
    except TypeError:
        save_info['\t L x L']: f'{param_set.L} x {param_set.L}'

    with open(f"{path_to_ensemble}/info/ensemble_info.txt", "a") as info_file:
        info_file.write("\n______Ensemble Parameters_______\n")
        for prm in save_info:
            info_file.write(prm + ' : ' + str(save_info[prm]) + '\n')
    return

def save_sim_out(path_to_ensemble:str, elapsed_time:str):
    """
    Save output info
    """
    with open(path_to_ensemble + "/info/ensemble_info.txt", "a") as info_file:
        info_file.write("\n______Out_______\n")
        info_file.write('\telapsed time : '+ elapsed_time + '\n')
        info_file.write("\tNotes : '...' ")  # Note section to document results
    return


def run_vel_ensemble(rho:float, beta:float, runs:int) -> dict:
    """
    Repeat velocity simulations over a number of `runs' for a given value of rho, beta.
    """
    import numpy as np
    ensemble_vel = np.zeros(runs)
    ensemble_perc = np.zeros(runs)
    settings = Settings()
    for N in range(runs):
        if settings.verbose>=2:
            print('Repeat : {}'.format(N))
        out_mts = model_dynamics.runSim(pc=ModelParamSet(rho, beta), metrics=Metrics(), settings=settings)[1]
        ensemble_vel[N] = out_mts.maxD.max()/out_mts.endT
        ensemble_perc[N] = perc(out_mts.percolation)
        if settings.verbose>=2:
            print("\t Percolation : {} @ t = {} (days)".format(out_mts.percolation, out_mts.endT))
            print('\t Max Distance : {} (m)'.format(round(out_mts.maxD.max(), 3)))
            print('\t ~ Spread Vel : {} (m/day)'.format(round(out_mts.maxD.max()/out_mts.endT, 3)))
            print('\t ~ Spread Vel: {} (km/ry)'.format(round((out_mts.maxD.max()*365)/(out_mts.endT*1000), 3)))
    return {'ensemble_velocity':ensemble_vel, 'ensemble_percolation':ensemble_perc}


def run_R0_ensemble(rho:float, beta:float, runs:int) -> dict:
    """
    Repeat R0 simulations over a number of `runs' for a given value of rho, beta.
    """
    import numpy as np
    from helper_methods import R0_generation_mean
    ensemble_R0 = []
    extinctionT = [None]*runs
    mortality_ratio = [None]*runs
    settings = Settings()
    for N in range(runs):
        if settings.verbose:
            print('Repeat : {}'.format(N))
        out_metrics = model_dynamics.runSim(pc=ModelParamSet(rho, beta), metrics=Metrics(), settings=settings)[1]
        R0_histories = R0_generation_mean(out_metrics.R0_histories)
        ensemble_R0.append(R0_histories)  # the number of gen-0 secondary infections
        extinctionT[N] = out_metrics.extinctionT
        mortality_ratio[N] = out_metrics.mortality_ratio
    return {'mean_R0_vs_gen_core_ensemble': ensemble_R0, 'extinction_time':extinctionT,
            'mortality_ratio':mortality_ratio}






