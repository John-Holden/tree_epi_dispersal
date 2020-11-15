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

def collect_data(name):
    vlist = sorted(os.listdir(name+'/vel/'))
    plist = sorted(os.listdir(name+'/vel/'))
    assert len(vlist) == len(plist)
    shape = np.load(name+'/vel/'+vlist[0])
    print(shape)
    for file in vlist:
        print(file)
    sys.exit()

ens_name = os.getcwd()+'/2020-11-14-hpc-vel-ens'
collect_data(ens_name)
sys.exit()
perc_ens_dat = np.load(ens_name+'/perc/perc_ensemble_.npy')
rhos = np.load(ens_name+'/info/rhos.npy')
betas = np.load(ens_name+'/info/betas.npy')
vel_ens_mean = vel_ens_dat.mean(axis=2)
perc_ens_mean = perc_ens_dat.mean(axis=2)
plot1D_mean(rhos, betas, vel_ens_mean, perc_ens_mean)
sys.exit()