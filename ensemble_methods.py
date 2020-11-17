import model
import sys
perc = lambda p: 1 if p else 0

def mk_new_dir(name):
    import os
    c=0
    for ith_name in sorted(os.listdir('./ensemble_dat/')):
        if name in ith_name:
            c+=1
    if c == 0:
        os.mkdir('./ensemble_dat/'+name)
        os.mkdir('./ensemble_dat/'+name+'/vel')
        os.mkdir('./ensemble_dat/'+name+'/perc')
        os.mkdir('./ensemble_dat/'+name+'/info')
        return name
    else:
        os.mkdir(os.getcwd()+'/ensemble_dat/'+name+'-{}'.format(c))
        os.mkdir('./ensemble_dat/'+name+'-{}/vel'.format(c))
        os.mkdir('./ensemble_dat/'+name+'-{}/perc'.format(c))
        os.mkdir('./ensemble_dat/'+name+'-{}/info'.format(c))
        return name+'-{}'.format(c)

def saveFunc(vel_ens_dat, perc_ens_dat, rhos, betas, mPs, sts, ensName, jobId):
    import numpy as np
    save_name = 'ensemble_dat/'+ensName
    np.save(save_name+"/vel/vel_ensemble_"+jobId, vel_ens_dat)
    np.save(save_name+"/perc/perc_ensemble_"+jobId, perc_ens_dat)
    if jobId == '1' or jobId == '':  # save ensemble details
        np.save(save_name+"/info/betas", betas)
        np.save(save_name+"/info/rhos", rhos)
        mPs = mPs(0, 0)
        save_info = {"model":mPs.model, "alpha": str(mPs.alpha) + '(m)',
                     "ell":str(mPs.ell)+'(m)', "L":str(mPs.L)+'X'+str(mPs.L),
                     "rhos":str([rhos[0], rhos[-1]])+' | '+str(len(rhos)),
                     "betas":str([betas[0], betas[-1]])+' | '+str(len(betas)),
                     "core repeats":vel_ens_dat.shape[2], "initial epi r": str(mPs.r),
                     "percolation boundary": sts.boundary}
        with open(save_name + "/info/ensemble_info.txt", "w+") as info_file:
            info_file.write("______Simulation Parameters_______" + "\n")
            for prm in save_info:
                info_file.write(prm + ' : ' + str(save_info[prm]) + '\n')
            info_file.write("Notes : '...' ")  # Note section to document results
    return

def runSIR_ensemble(N, mPrm, mts, sts) -> 'Success':
    import pickle
    from plots.plotLib import pltSIR
    ensemble_results = {}
    for i in range(N):
        print('Repeat : {}'.format(i))
        [out_mPrm, out_mts] = model.runSim(pc=mPrm(), metrics=mts(), settings=sts())
        ensemble_results["S{}".format(i)] = out_mts.numS
        ensemble_results["I{}".format(i)] = out_mts.numI
        ensemble_results["R{}".format(i)] = out_mts.numR
        pltSIR(S=out_mts.numS, I=out_mts.numI, R=out_mts.numR, dt=1)

    file = open("./SIR/ENSEMBLE.dat", "wb")
    pickle.dump(ensemble_results, file)
    file.close()
    return "Success"

def runVel_ensemble(r, b, runs, MPrm, Mts, Sts) -> '[velocity ensemble, percolation ensemble]':
    """
    :param r: float, rho av tree density
    :param b: float, beta infectivity constant
    :param N: int, number of repeats
    :param mPrm: class, model parameters
    :param mts: class, metrics & time-series
    :param sts: class, model settings
    :return: arr (1, N), array of repeats
    """
    import numpy as np
    ensemble_vel = np.zeros(runs)
    ensemble_perc = np.zeros(runs)
    sts = Sts()
    for N in range(runs):
        if sts.verbose2:
            print('Repeat : {}'.format(N))
        [out_mPrm, out_mts] = model.runSim(pc=MPrm(r, b), metrics=Mts(), settings=sts)
        ensemble_vel[N] = out_mts.maxD.max()/out_mts.endT
        ensemble_perc[N] = perc(out_mts.percolation)
        if sts.verbose2:
            print("\t Percolation : {} @ t = {} (days)".format(out_mts.percolation, out_mts.endT))
            print('\t Max Distance : {} (m)'.format(round(out_mts.maxD.max(), 3)))
            print('\t ~ Spread Vel : {} (m/day)'.format(round(out_mts.maxD.max()/out_mts.endT, 3)))
            print('\t ~ Spread Vel: {} (km/ry)'.format(round((out_mts.maxD.max()*365)/(out_mts.endT*1000), 3)))
    return ensemble_vel, ensemble_perc



