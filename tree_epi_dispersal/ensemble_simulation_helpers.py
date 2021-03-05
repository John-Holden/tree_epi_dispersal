import os
import math
import numpy as np
from typing import Union
from warnings import warn

from parameters_and_settings import ModelParamSet, Settings, Metrics

perc = lambda p: 1 if p else 0

def time_print(time_seconds:int, msg:str='Simulation done in: ') -> str:
    """
    Pretty formatter to display time
    """
    seconds = math.floor(time_seconds)
    hrs = math.floor(seconds / 3600)
    mns = math.floor(seconds / 60)
    secs = seconds % 60

    if seconds < 60:
        print(f'{msg} {seconds} (s)')
    elif 60 <= seconds < 3600:
        print(f'{msg} {mns} (mins): {seconds} (s)')
    elif seconds >= 3600:
        print(f'{msg} {hrs} (Hrs): {mns%60} (mins): {secs} (s)')

    return f'{msg} {hrs} (Hrs): {mns%60} (mins): {secs} (s)'


def mk_new_dir(name: str):
    "Save new directory to file."

    dir1 = f'{os.getcwd()}/temp_dat_store/{name}'
    sub_dir1 = f'{os.getcwd()}/temp_dat_store/{name}/info/'
    sub_dir2 = f'{os.getcwd()}/temp_dat_store/{name}/core_output'

    if os.path.exists(dir1):
        msg = f'FileExistsWarn: {os.getcwd()}/temp_dat_store/{name}'
        # raise FileExistsError(F'{os.getcwd()}/temp_dat_store/{name}')
        warn(msg)
        assert os.path.exists(sub_dir1)
        assert os.path.exists(sub_dir2)
        return

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

def save_meta_data(path_to_ensemble:str, elapsed_time:str):
    """
    Save output info
    """

    print(vars(ModelParamSet(0, 0)))

    assert 0
    with open(path_to_ensemble + "/info/ensemble_info.txt", "a") as info_file:
        info_file.write("\n______Out_______\n")
        info_file.write('\telapsed time : '+ elapsed_time + '\n')
        info_file.write("\tNotes : '...' ")  # Note section to document results
    return







