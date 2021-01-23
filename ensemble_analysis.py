"""
Plot and analyse ensemble data.
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt


pltParams = {'figure.figsize': (7.5, 5.5),
             'axes.labelsize': 15,
             'ytick.labelsize': 20,
             'xtick.labelsize': 20,
             'legend.fontsize': 'x-large'}
plt.rcParams.update(pltParams)

def plot1D_mean(rhos, betas, ens_mean, save=False):
    for i, beta in enumerate(betas):
        plt.plot(rhos, ens_mean[i], c='C'+str(i),
                 label=r'$\beta$ = {}'.format(round(beta, 8)))
        plt.scatter(rhos, ens_mean[i], c='C'+str(i),)

    plt.hlines(y=1, xmin=rhos[0], xmax=rhos[-1], ls='--')
    plt.legend()
    if save:
        plt.savefig('ens_dat.pdf')
    plt.show()
    return "SUCCESS"

def collect_data(name: str, field: str) -> 'np.ndarray | ensemble average of field f':
    """
    Collect each core result, as a np.ndarray, and average.
    """
    # todo test this and streamline
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
    plot1D_mean(rhos, betas, ens_mean, save=True)
    return ens_mean


def process_R0_ensemble(path_to_ensemble:str, ensemble_mean:dict):
    from helper_methods import avg_multi_dim
    box_sizes = np.load(f'{path_to_ensemble}/info/box_sizes.npy')

    for i, box_size in enumerate(ensemble_mean):
        assert int(box_size) == box_sizes[i]
        plt.plot(avg_multi_dim(ensemble_mean[box_size]), label=f' size: {box_size}')

    plt.legend()
    plt.show()


def ens_avg_dict_of_arrays(path_to_ensemble:str, metric:str) -> dict:
    import json
    from collections import defaultdict
    from helper_methods import avg_multi_dim
    """
    Iteratively load json object from a directory, each object represents a core. Average each core and return.
    """
    f_list = sorted(os.listdir(f'{path_to_ensemble}/{metric}/'))
    core_means = defaultdict(list)
    for core_result_name in f_list:
        with open(f"{path_to_ensemble}/{metric}/{core_result_name}") as f:
            core_R0_history = json.load(f)
            for box_size in core_R0_history:
                core_means[box_size].append(avg_multi_dim(core_R0_history[box_size]))
    return core_means



ens_name = os.getcwd()+'/data_store/2021-01-23-local-ensemble'
ensemble_avg = ens_avg_dict_of_arrays(ens_name, metric='R0_histories')
process_R0_ensemble(ens_name, ensemble_avg)