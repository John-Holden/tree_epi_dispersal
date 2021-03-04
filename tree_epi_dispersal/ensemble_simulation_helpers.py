import os
import numpy as np
from typing import Union
from tree_epi_dispersal import model_dynamics
from parameters_and_settings import Metrics, ModelParamSet, Settings

perc = lambda p: 1 if p else 0

def mk_new_dir(name: str):
    "Save new directory to file."

    if os.path.exists(f'{os.getcwd()}/temp_dat_store/{name}'):
        raise FileExistsError(f'{os.getcwd()}/temp_dat_store/{name}')

    os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}')
    os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/info/')
    os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/core_output')


def save_ens_info(ens_field_names: list, rhos: Union[np.ndarray, float], betas:Union[np.ndarray, float],
                  path_to_ensemble:str, per_core_repeats:int, box_sizes=None):
    """
    Save simulation ensemble to file at the beginning of the simulation.
    """
    from datetime import datetime
    param_set = ModelParamSet(0, 0)  # init
    save_info = {"\t ensemble averaged metrics: ": ens_field_names,
                 "\t start time: ": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                 "\t dispersal model: ": param_set.model,
                 "\t scale parameter $alpha$: ": f'{param_set.alpha} (m)',
                 "\t dispersal constant $ell$: ": f'{param_set.ell} (m)',
                 "\t infectious life-time $T$": f'{param_set.infLT}',
                 "\t $I -> R$ transitions: ": ' Exp ',
                 "\t # core repeats: ": f'{per_core_repeats}',
                 "\t initial epicenter radius: ": f'{param_set.r}',
                 "\t percolation boundary: ": f'{Settings.boundary}'}

    try:  # save iterable parameter types.
        iter(rhos)
        np.save(f"{path_to_ensemble}/info/rhos", rhos)
        save_info['\t rhos'] = f'{rhos} | {len(rhos)}'
    except TypeError:
        save_info["\t rho"] = f'{round(rhos, 5)}'
    try:
        iter(betas)
        np.save(f"{path_to_ensemble}/info/betas", betas)
        save_info["\t betas"] = f'{betas} | {len(betas)}'
    except TypeError:
        save_info["\t beta"] = f'{round(betas, 7)}'
    try:
        iter(box_sizes)
        np.save(f"{path_to_ensemble}/info/box_sizes", box_sizes)
        save_info["\t domain size "] = f'{box_sizes}'
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
    from tree_epi_dispersal.model_dynamics_helpers import R0_generation_mean
    ensemble_R0 = []
    for N in range(runs):
        if Settings.verbose:
            print('Repeat : {}'.format(N))
        out_metrics = model_dynamics.runSim(pc=ModelParamSet(rho, beta), metrics=Metrics())[1]
        R0_histories = R0_generation_mean(out_metrics.R0_histories)
        ensemble_R0.append(R0_histories)  # the number of gen-0 secondary infections
    return {'mean_R0_vs_gen_core_ensemble': ensemble_R0}






