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

def ens_avg_dict_of_arrays(name:str, metric:str):
    """

    """
    import pickle
    f_list = sorted(os.listdir(f'{name}/{metric}/'))
    for core_result_name in f_list:
        core_ens = open(f"{name}/{metric}/{core_result_name}", 'rb')
        for key in core_ens:
            print(key)
        sys.exit()
    # object_file = pickle.load(file)


    # for box_size in R0_gen_ens:
    #     plt.plot(avg_multi_dim(R0_gen_ens[box_size]), label=f' size: {box_sizes[box_size]}, '
    #                                                         f'grid size:{box_sizes[box_size]*mprms.alpha/1000}km')
    # plt.xlabel('Generation', size=20)
    # plt.ylabel(r'$\overline{R}_0$', size=25)
    # plt.legend()
    # plt.show()



ens_name = os.getcwd()+'/data_store/2021-01-23-hpc-R0-generation'
ens_avg_dict_of_arrays(ens_name, metric='R0_histories')