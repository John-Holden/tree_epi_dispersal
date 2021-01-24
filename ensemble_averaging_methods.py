import model
from typing import Union, Type
import numpy as np
perc = lambda p: 1 if p else 0

def mk_new_dir(name):
    import os
    import sys
    if os.path.exists(f'{os.getcwd()}/ensemble_dat/{name}'):
        sys.exit(f'ERROR DUPLICATE DIRECTORY{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/vel')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/R0_histories')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/extinction_time')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/mortality_ratio')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/perc')
    os.mkdir(f'{os.getcwd()}/ensemble_dat/{name}/info')
    return name


def save_ensemble(ens_results, ens_field_names, ensemble_name, job_id) -> 'Success':
    """
    Save the ensemble results
    :param ens_results: list [arg1,...], arg1: npy.arr
    :param ens_field_names: list [name1,..] name1:str
    :param ensemble_name: str:
    :param job_id: str:
    :return: Success
    """
    save_name = 'ensemble_dat/'+ensemble_name
    assert len(ens_field_names) == len(ens_results), "/error. supplied field names must be same length" \
                                                     " as saved ensemble filed."
    if job_id == '':  # local mode
        job_id='lcl_test'
    for c, metric in enumerate(ens_field_names):
        np.save(save_name + f"/{metric}/" + job_id, ens_results[c])
    return 'Success'

def save_ens_info(ens_field_names: list, rhos: Union[np.ndarray, float], betas:Union[np.ndarray, float],
                  param_set, settings, ensemble_name, per_core_repeats:int, box_sizes:Union[None, list]):
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
                 "\t L": str(param_set.L) + ' X ' + str(param_set.L),
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


def runSIR_ensemble(N, mPrm, mts, sts) -> 'Success':
    import pickle
    from plots.plotLib import pltSIR
    ensemble_results = {}
    for i in range(N):
        print('Repeat : {}'.format(i))
        [out_mPrm, out_mts] = model.runSim(pc=mPrm(), metrics=mts(), settings=sts())
        ensemble_results["S{}".format(i)] = out_mts.numS
        ensemble_results["I{}".format(i)] = out_mts.numI
        ensemble_results["R{}".format(i)] = out_mts.numR
        pltSIR(S=out_mts.numS, I=out_mts.numI, R=out_mts.numR, dt=1)

    file = open("./SIR/ENSEMBLE.dat", "wb")
    pickle.dump(ensemble_results, file)
    file.close()
    return "Success"


def runVel_ensemble(r, b, runs, MPrm, Mts, Sts) -> '[velocity ensemble, percolation ensemble]':
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
        [out_mPrm, out_mts] = model.runSim(pc=MPrm(r, b), metrics=Mts(), settings=sts)
        ensemble_vel[N] = out_mts.maxD.max()/out_mts.endT
        ensemble_perc[N] = perc(out_mts.percolation)
        if sts.verbose:
            print("\t Percolation : {} @ t = {} (days)".format(out_mts.percolation, out_mts.endT))
            print('\t Max Distance : {} (m)'.format(round(out_mts.maxD.max(), 3)))
            print('\t ~ Spread Vel : {} (m/day)'.format(round(out_mts.maxD.max()/out_mts.endT, 3)))
            print('\t ~ Spread Vel: {} (km/ry)'.format(round((out_mts.maxD.max()*365)/(out_mts.endT*1000), 3)))
    return ensemble_vel, ensemble_perc


def runR0_ensemble(r, b, runs, MPrm, Mts, Sts) -> '[R0, extinctinon time, mortality ratio]':
    """
    :param r: float, rho
    :param b: float, beta
    :param runs: int
    :param MPrm: Class, model parameters
    :param Mts: Class, metrics
    :param Sts: Class, settings
    """
    import numpy as np
    from helper_methods import R0_generation_mean
    ensemble_R0 = np.zeros(runs)
    extinctionT = np.zeros(runs)
    mortality_ratio = np.zeros(runs)
    sts = Sts()
    for N in range(runs):
        if sts.verbose:
            print('Repeat : {}'.format(N))
        [out_mPrm, out_mts] = model.runSim(pc=MPrm(r, b), metrics=Mts(), settings=sts)
        if len(R0_generation_mean(out_mts.R0_trace)) == 0:
            ensemble_R0[N] = 0
        else:
            ensemble_R0[N] = R0_generation_mean(out_mts.R0_trace).mean()
        extinctionT[N] = out_mts.extinctionT
        mortality_ratio[N] = out_mts.mortality_ratio

    return [ensemble_R0, extinctionT, mortality_ratio]






