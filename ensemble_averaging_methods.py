import model_dynamics
from typing import Union, Type
import numpy as np
perc = lambda p: 1 if p else 0
from PARAMETERS_AND_SETUP import Metrics, ModelParamSet, Settings

def mk_new_dir(name):
    import os
    import sys
    if os.path.exists(f'{os.getcwd()}/ensemble_dat/{name}'):
        sys.exit(f'ERROR DUPLICATE DIRECTORY{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/info/')
    return name


def save_ens_info(ens_field_names: list, rhos: Union[np.ndarray, float], betas:Union[np.ndarray, float],
                  param_set, settings, ensemble_name, per_core_repeats:int, box_sizes=None):
    """
    Save simulation ensemble to file at the beginning of the simulation.
    """
    from datetime import datetime
    save_name = 'ensemble_dat/'+ensemble_name
    save_info = {"\t ensemble averaged metrics: ": ens_field_names,
                 "\t start time: ": datetime.now(),
                 "\t model: ": param_set.model,
                 "\t alpha: ": str(param_set.alpha) + '(m)',
                 "\t ell: ": str(param_set.ell) + '(m)',
                 "\t core repeats: ": per_core_repeats,
                 "\t initial epicenter radius: ": str(param_set.r),
                 "\t percolation boundary: ": settings.boundary}

    # save iterable parameter types.
    if type(rhos) is np.ndarray:
        np.save(save_name+"/info/betas", betas)
        save_info['\t rhos'] = f'{rhos[0]}, {rhos[-1]} | {len(rhos)}'
    elif type(rhos) is float:
        save_info["\t rho"] = f'{round(rhos, 5)}'
    if type(betas) is np.ndarray:
        np.save(save_name+"/info/betas", betas)
        save_info["\t betas"] = f'{betas[0]}, {betas[-1]} | {len(betas)}'
    elif type(betas) is float:
        save_info["\t beta"] = f'{round(betas, 7)}'
    if type(box_sizes) is list or type(box_sizes) is np.ndarray:
        np.save(save_name+"/info/box_sizes", box_sizes)
        save_info["\t domain size "] = [L for L in box_sizes]
        save_info["\t modelled size"] = [f'{L*param_set.alpha/1000} x {L*param_set.alpha/1000} (km)' for L in box_sizes]
    elif box_sizes is None:
        save_info['\t L x L']: f'{param_set.L} x {param_set.L}'

    with open(save_name + "/info/ensemble_info.txt", "a") as info_file:
        info_file.write("\n______Ensemble Parameters_______\n")
        for prm in save_info:
            info_file.write(prm + ' : ' + str(save_info[prm]) + '\n')
    return

def save_sim_out(ensemble_name:str, elapsed_time:str):
    """
    Save output info
    """
    save_name = 'ensemble_dat/' + ensemble_name
    with open(save_name + "/info/ensemble_info.txt", "a") as info_file:
        info_file.write("\n______Out_______\n")
        info_file.write('\telapsed time : '+ elapsed_time + '\n')
        info_file.write("\tNotes : '...' ")  # Note section to document results
    return


def run_vel_ensemble(r, b, runs, MPrm, Mts, Sts) -> '[velocity ensemble, percolation ensemble]':
    """
    :param r: float, rho av tree density
    :param b: float, beta infectivity constant
    :param N: int, number of repeats
    :param mPrm: class, model parameters
    :param mts: class, metrics & time-series
    :param sts: class, model settings
    :return: arr (1, N), array of repeats
    """
    import numpy as np
    ensemble_vel = np.zeros(runs)
    ensemble_perc = np.zeros(runs)
    sts = Sts()
    for N in range(runs):
        if sts.verbose:
            print('Repeat : {}'.format(N))
        [out_mPrm, out_mts] = model_dynamics.runSim(pc=MPrm(r, b), metrics=Mts(), settings=sts)
        ensemble_vel[N] = out_mts.maxD.max()/out_mts.endT
        ensemble_perc[N] = perc(out_mts.percolation)
        if sts.verbose:
            print("\t Percolation : {} @ t = {} (days)".format(out_mts.percolation, out_mts.endT))
            print('\t Max Distance : {} (m)'.format(round(out_mts.maxD.max(), 3)))
            print('\t ~ Spread Vel : {} (m/day)'.format(round(out_mts.maxD.max()/out_mts.endT, 3)))
            print('\t ~ Spread Vel: {} (km/ry)'.format(round((out_mts.maxD.max()*365)/(out_mts.endT*1000), 3)))
    return ensemble_vel, ensemble_perc


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






