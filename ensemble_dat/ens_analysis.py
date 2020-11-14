import sys
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

# ens_name = './2020-11-14-hpc-vel-ens'
ens_name = './2020-11-14-lcl-vel-ens-1'
vel_ens_dat = np.load(ens_name+'/vel/vel_ensemble_.npy')
perc_ens_dat = np.load(ens_name+'/perc/perc_ensemble_.npy')
rhos = np.load(ens_name+'/info/rhos.npy')
betas = np.load(ens_name+'/info/betas.npy')
vel_ens_mean = vel_ens_dat.mean(axis=2)
perc_ens_mean = perc_ens_dat.mean(axis=2)
plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean)
sys.exit()