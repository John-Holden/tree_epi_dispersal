import sys, os
import numpy as np
import matplotlib.pyplot as plt

def plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean):
    for i, beta in enumerate(betas):
        plt.plot(rhos, vel_ens_mean[i], c='C'+str(i),
                 label='beta = {}'.format(round(beta, 3)))
        # plt.plot(rhos, vel_ens_mean[i] * perc_ens_mean[i],
        #          c='C'+str(i), alpha=0.50, ls='--')
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

ens_name = os.getcwd()+'/data_store/2020-11-15-hpc-vel-ens'
rhos = np.load(ens_name+'/info/rhos.npy')
betas = np.load(ens_name+'/info/betas.npy')
vel_ens_mean = collect_data(ens_name, 'vel')
perc_ens_mean = collect_data(ens_name, 'perc')
plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean)
sys.exit()