import sys, os
import numpy as np
import matplotlib.pyplot as plt

def plot1D_mean(rhos, betas, ens_mean):
    for i, beta in enumerate(betas):
        plt.plot(rhos, ens_mean[i], c='C'+str(i),
                 label='beta = {}'.format(round(beta, 8)))
        plt.scatter(rhos, ens_mean[i], c='C'+str(i),)

    plt.hlines(y=1, xmin=rhos[0], xmax=rhos[-1], ls='--')
    plt.legend()
    plt.show()
    return "SUCCESS"

def collect_data(name, field) -> 'ensemble average of field f':
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

def run_plot(name, metric):
    rhos = np.load(name + '/info/rhos.npy')
    betas = np.load(name + '/info/betas.npy')
    ens_mean = collect_data(name, metric)
    np.save('R0_data', ens_mean)
    plot1D_mean(rhos, betas, ens_mean)
    return

ens_name = os.getcwd()+'/data_store/2020-11-26-hpc-R0-trace'
run_plot(ens_name, metric='R0_trace')
sys.exit()