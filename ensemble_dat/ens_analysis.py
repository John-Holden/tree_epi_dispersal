import sys, os
import numpy as np
import matplotlib.pyplot as plt

def plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean):
    for i, beta in enumerate(betas):
        print('beta : ', vel_ens_mean[i])
        plt.plot(rhos, vel_ens_mean[i],
                 label='beta = {}'.format(round(beta, 3)))
        plt.plot(rhos, vel_ens_mean[i] * perc_ens_mean[i],
                 alpha=0.50, ls='--')
    plt.legend()
    plt.show()
    return "SUCCESS"

def collect_data(name, field) -> 'ensemble average of field f':
    f_list = sorted(os.listdir(name+'/'+field+'/'))
    dat = np.load(name+'/'+field+'/'+f_list[0])
    print('\t @collecting data: field {}'.format(field))
    for file in f_list[1:]:
        print(file)
        dat += np.load(name+'/'+field+'/'+file)
    dat = dat / len(f_list)
    return dat.mean(axis=2)

ens_name = os.getcwd()+'/2020-11-15-hpc-vel-ens'
rhos = np.load(ens_name+'/info/rhos.npy')
betas = np.load(ens_name+'/info/betas.npy')
vel_ens_mean = collect_data(ens_name, 'vel')
perc_ens_mean = collect_data(ens_name, 'perc')

plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean)
sys.exit()