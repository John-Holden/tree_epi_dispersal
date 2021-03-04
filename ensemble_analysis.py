"""
Plot and analyse ensemble data.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

PATH_TO_DATA_STORE = os.getcwd()+'/data_store'

pltParams = {'figure.figsize': (7.5, 5.5),
             'axes.labelsize': 15,
             'ytick.labelsize': 20,
             'xtick.labelsize': 20,
             'legend.fontsize': 'x-large'}
plt.rcParams.update(pltParams)


def plot_rho_beta_ensemble_1D(ensemble:np.ndarray, rhos:np.ndarray, betas:np.ndarray):
    """
    For an ensemble, plot a series of 1D lines, rho axis, for multiple beta values.
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    for i in range(len(ensemble)):
        rho_line = ensemble[i]
        ax.plot(rhos, rho_line, label=f'beta = {round(betas[i], 5)}')
        ax.scatter(rhos, rho_line)
    print('BETAS: ', betas)
    ax.plot([0, rhos[-1]], [1, 1], color='r', ls='--')

    plt.legend()
    plt.show()


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


def plot_R0_ens_vs_gen(path_to_ensemble:str, ensemble_mean:dict, save=False):
    """
    Process ensemble_core_averaged R0 history and plot generational mean R0
    """
    from helper_methods import avg_multi_dim
    box_sizes = np.load(f'{path_to_ensemble}/info/box_sizes.npy')[::-1]
    c = 0
    for box_size in ensemble_mean:
        assert int(box_size) in box_sizes
        c += 1
    assert c == len(box_sizes)
    R0_v_L = np.zeros(len(box_sizes))

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for i, box_size in enumerate(box_sizes):
        R0_vs_gen = avg_multi_dim(ensemble_mean[str(box_size)])[:21]
        R0_v_L[i] = R0_vs_gen[0]
        gen = np.arange(1, len(R0_vs_gen)+1)
        ax.plot(gen, R0_vs_gen, label=f'{box_size}  {box_size}')
        ax.scatter(gen, R0_vs_gen)

    ax.set_xlim(0.5, 10.5)
    plt.legend()
    if save:
        np.save('R0_vs_L_alpha_()m', R0_v_L)
        plt.tight_layout()
        plt.savefig('R0_vs_gen.pdf')
    plt.show()


def plot_R0_ens_vs_L(save=False):
    """
    For different alpha values, plot saved R0 values against domain size L
    """
    data_sets = {'2021-01-24-hpc-R0-generation-alpha-5': '5',
                 '2021-01-24-hpc-R0-generation-alpha-10': '10',
                 '2021-01-24-hpc-R0-generation-alpha-7_5': '7_5'}

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for data_set in data_sets:
        box_sizes = np.load(f'{PATH_TO_DATA_STORE}/{data_set}/info/box_sizes.npy')[::-1]
        alpha_5 = np.load(f'{PATH_TO_DATA_STORE}/{data_set}/R0_vs_L_alpha_{data_sets[data_set]}m.npy')
        ax.plot(box_sizes, alpha_5, label=rf'$\alpha = $ {data_sets[data_set].replace("_", ".")}')
        ax.scatter(box_sizes, alpha_5)

    plt.legend()
    plt.tight_layout()
    if save:
        plt.savefig('R0_vs_L_vs_alpha.pdf')
    plt.show()


def ens_avg_dict_of_R0_arrays(path_to_ensemble:str, metric:str) -> dict:
    """
    Iteratively load json object from a directory, each object represents a core. Average each core and return.
    """
    import json
    from collections import defaultdict
    from helper_methods import avg_multi_dim

    f_list = sorted(os.listdir(f'{path_to_ensemble}/{metric}/'))
    core_means = defaultdict(list)
    for core_result_name in f_list:
        with open(f"{path_to_ensemble}/{metric}/{core_result_name}") as f:
            core_R0_history = json.load(f)
            for box_size in core_R0_history:
                core_means[box_size].append(avg_multi_dim(core_R0_history[box_size]))


    print(f'Ensemble size {len(f_list) * len(core_means[box_size])}')
    return core_means


def process_avg_R0(R0_struct:list) -> float:
    R0_av = 0
    for R0_vs_gen in R0_struct:
        if len(R0_vs_gen):  # zero-length
            R0_av += R0_vs_gen[0]  # sum first-generation R0

    return R0_av/len(R0_struct)

def write_package(ensemble: np.ndarray, path_to_ens:str):
    """
    Create folder to be used for land-scape control code-base.
    """
    path_to_save = f'{path_to_ens}/landscape_control_package'
    if os.path.exists(f'{path_to_ens}/landscape_control_input'):
        print(f'Warning, folder {path_to_ens}/landscape_control_input already exists!')
        return
    else:
        import shutil
        os.mkdir(path_to_save)
        np.save(f'{path_to_save}/ensemble', ensemble)
        shutil.copy(f'{path_to_ens}/info/ensemble_info.txt', f'{path_to_save}/ensemble_info.txt')
        shutil.copy(f'{path_to_ens}/info/rhos.npy', f'{path_to_save}/rhos.npy')
        shutil.copy(f'{path_to_ens}/info/betas.npy', f'{path_to_save}/betas.npy')


def process_R0_ensemble(path_to_ensemble:str, field_of_interest:str,
                        produce_landscape_control_package:bool) -> np.ndarray:
    """Load json, for each rho-beta key find R0, then average over of all core results. """
    import json
    rhos = np.load(f'{path_to_ensemble}/info/rhos.npy')
    betas = np.load(f'{path_to_ensemble}/info/betas.npy')

    if os.path.exists(f'{path_to_ensemble}/R0-vs-rho.npy'): # File already processed.
        ensemble = np.load(f'{path_to_ensemble}/R0-vs-rho.npy')

    else:  # Iterate through all core output, collect and average.
        f_list = sorted(os.listdir(f'{path_to_ensemble}/core_output/'))
        with open(f'{path_to_ensemble}/info/ensemble_info.txt') as ens_info:
            for line in ens_info.readlines():
                if 'core repeats' in line:
                    number_of_core_repeats = int(line.split()[-1])
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
                        R0_vs_gen_v_ens = core_result[f'rho_{rho}_beta_{beta}'][field_of_interest]   # list of lists
                        ensemble[i, j] += process_avg_R0(R0_vs_gen_v_ens)  # find avg R0 for (rho, beta)

        ensemble = ensemble/len(f_list)

    if produce_landscape_control_package:
        write_package(ensemble, path_to_ensemble)

    return ensemble, rhos, betas

if __name__ == '__main__':
    ens_name = f'{PATH_TO_DATA_STORE}/2021-01-30-hpc-R0-vs-rho'
    ensemble, rhos, betas = process_R0_ensemble(ens_name, field_of_interest='mean_R0_vs_gen_core_ensemble',
                                                produce_landscape_control_package=True)

    plot_rho_beta_ensemble_1D(ensemble, rhos, betas)
