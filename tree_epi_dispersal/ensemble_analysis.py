import os
import json
import numpy as np
from typing import Union


def collect_data(name: str, field: str) -> 'np.ndarray | ensemble average of field f':
    """
    Collect each core result, as a np.ndarray, and average.
    """
    f_list = sorted(os.listdir(name+'/'+field+'/'))
    dat = np.load(name+'/'+field+'/'+f_list[0])
    print('@collecting data: field {}'.format(field))
    for file in f_list[1:]:
        dat += np.load(name+'/'+field+'/'+file)

    dat = dat / len(f_list)
    print('\t cores = {} '.format(len(f_list)))
    print('\t repeats/cores = {} '.format(dat.shape[2]))
    print('\t -> ensemble size = {} '.format(dat.shape[2] * len(f_list)))
    return dat.mean(axis=2)


def collect_and_plot(name: str, metric: str):
    """
    Given the dataset name and the metric of interest, load np.array-based ensemble data (n-dimensional)
    for rhos, betas and save average in ensemble directory.
    """
    rhos = np.load(name + '/info/rhos.npy')
    betas = np.load(name + '/info/betas.npy')
    ens_mean = collect_data(name, metric)
    np.save(name + '/ens_R0_data', ens_mean)
    return ens_mean, rhos, betas


def ens_avg_dict_of_R0_arrays(path_to_ensemble:str, metric:str) -> dict:
    """
    Iteratively load json object from a directory, each object represents a core. Average each core and return.
    """
    import json
    from collections import defaultdict
    from tree_epi_dispersal.model_dynamics_helpers import avg_multi_dim

    f_list = sorted(os.listdir(f'{path_to_ensemble}/{metric}/'))
    core_means = defaultdict(list)
    for core_result_name in f_list:
        with open(f"{path_to_ensemble}/{metric}/{core_result_name}") as f:
            core_R0_history = json.load(f)
            for box_size in core_R0_history:
                core_means[box_size].append(avg_multi_dim(core_R0_history[box_size]))

    print(f'Ensemble size {len(f_list) * len(core_means[box_size])}')
    return core_means


def process_avg_R0_struct(R0_struct:dict, gen:Union[None, int] = None):
    """"From the R0-history dictionary, process statistics of each infectious tree."""
    R0_cumulative = np.zeros(1000)
    counts = np.zeros(1000)
    max_gen = 0
    for site, R0_statistics in R0_struct.items():
        R0_cumulative[R0_statistics[1]] += R0_statistics[0]
        counts[R0_statistics[1]] += 1
        max_gen = R0_statistics[1] if R0_statistics[1] > max_gen else max_gen

    R0_cumulative = R0_cumulative[:max_gen + 1]
    counts = counts[:max_gen + 1]

    if gen is None:
        return R0_cumulative / counts

    return R0_cumulative[gen-1] / counts[gen-1]


def R0_statistics_selector(R0_gens: list, nth_gen: int = 1):
    """
    Return the desired statistics from the R0_vs_gen data.
    :param R0_gens: [ avg(R0_1), avg(R0_2),...,avg(R0_N)]
    :param nth_gen: the type of stats we desire from the processed ensemble
    :return:
    """
    return sum(R0_gens[:nth_gen]) / len(R0_gens[:nth_gen])


def process_avg_R0(R0_struct: list, process_nth_gen: int) -> float:
    R0_av = 0
    for R0_vs_gen in R0_struct:
        if len(R0_vs_gen):  # zero-length
            R0_av += R0_statistics_selector(R0_vs_gen, process_nth_gen)  # sum upto n^th-generation R0

    return R0_av/len(R0_struct)


def write_package(ensemble: np.ndarray, path_to_ens: str):
    """
    Create folder to be used for land-scape control code-base.
    """
    path_to_save = f'{path_to_ens}/landscape_control_package'
    if os.path.exists(f'{path_to_ens}/landscape_control_package'):
        print(f'Warning, folder {path_to_ens}/landscape_control_input already exists!')
        return
    else:
        import shutil
        os.mkdir(path_to_save)
        np.save(f'{path_to_save}/ensemble', ensemble)
        shutil.copy(f'{path_to_ens}/info/ensemble_info.txt', f'{path_to_save}/ensemble_info.txt')
        shutil.copy(f'{path_to_ens}/info/rhos.npy', f'{path_to_save}/rhos.npy')
        shutil.copy(f'{path_to_ens}/info/betas.npy', f'{path_to_save}/betas.npy')


def process_R0_ensemble(path_to_ensemble: str, produce_landscape_control_package: bool,
                        process_nth_gen: int = 1) -> np.ndarray:
    """Load json, for each rho-beta key find R0, then average over of all core results. """
    rhos = np.load(f'{path_to_ensemble}/info/rhos.npy')
    betas = np.load(f'{path_to_ensemble}/info/betas.npy')

    if os.path.exists(f'{path_to_ensemble}/R0-vs-rho.npy'):  # File already processed.
        ensemble = np.load(f'{path_to_ensemble}/R0-vs-rho.npy')

    else:  # Iterate through all core output, collect and average.
        f_list = sorted(os.listdir(f'{path_to_ensemble}/core_output/'))
        with open(f'{path_to_ensemble}/info/ensemble_info.txt') as ens_info:
            for line in ens_info.readlines():
                if 'ensemble_size' in line:
                    number_of_core_repeats = int(line.split(' ')[-1])

        ensemble_size = len(f_list) * number_of_core_repeats

        print(f'Ensemble size = {ensemble_size}')

        with open(f'{path_to_ensemble}/info/ensemble_info.txt', 'a') as ens_info:  # write total size to file
            ens_info.write(f'Ensemble size : {ensemble_size}')

        ensemble = np.zeros(shape=[ len(betas), len(rhos)])
        for core_result_name in f_list:
            with open(f"{path_to_ensemble}/core_output/{core_result_name}") as f:
                core_result = json.load(f)
                for i, beta in enumerate(betas):
                    for j, rho in enumerate(rhos):
                        R0_vs_gen_v_ens = core_result[f'rho_{rho}_beta_{beta}']
                        ensemble[i, j] += process_avg_R0(R0_vs_gen_v_ens, process_nth_gen)  # find avg R0 for (rho, beta)

        ensemble = ensemble/len(f_list)

    if produce_landscape_control_package:

        write_package(ensemble, path_to_ensemble)

    return ensemble, rhos, betas