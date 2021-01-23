import numpy as np
from typing import Type
from helper_methods import set_R0trace, ijDistance, get_new_I, setFields
from run_simulation import ModelParamSet, Settings, Metrics

printStep = lambda t, freq : print('\t\t Time : {} (days)'.format(t)) if t % freq == 0 else None

def runSim(pc: Type[ModelParamSet], metrics: Type[Metrics], settings: Type[Settings]) -> '[S,I,R], metrics':
    """
    Run dynamic simulation of pathogen spread,
    :return: SIR fields and updated metrics
    """
    S,I,R = setFields(pc.L, pc.rho, pc.epiC, pc.r, pc.init_n_infected)

    if settings.plot:
        from plots.plotLib import pltSim

    for t in range(pc.tend):
        if settings.verbose == 2:
            printStep(t, freq=1)
        S_ = np.where(S)
        I_ = np.where(I)
        R_ = np.where(R)
        # update metrics
        if t == 0:
            set_R0trace(I_, metrics.R0_histories)
        metrics.numS[t] = len(S_[0])
        metrics.numI[t] = len(I_[0])
        metrics.numR[t] = len(R_[0])
        if metrics.numI[t] == 0:  # BCD 1
            metrics.extinction = True
            metrics.extinctionT = t
            break

        metrics.maxD[t] = ijDistance(i=[pc.epiC, pc.epiC], j=I_).max() * pc.alpha
        if not metrics.percolation:
            if metrics.maxD.max() >= (pc.L/2 - 10) * pc.alpha:
                metrics.percT = t
                metrics.percolation = True

        if settings.boundary and metrics.percolation:
            break
        # update fields S, I, R
        newI_ind, max_gen_exceeded = get_new_I(S_, I_, pc.beta, pc.alpha, pc.ell, pc.model,
                                                   metrics.R0_histories, settings.gen_limit)

        if settings.gen_limit is not None and max_gen_exceeded:
            # if no remaining infected trees of order `gen-limit', terminate simulation
            break

        S[newI_ind] = 0
        I[newI_ind] = 1
        I = I + (I>=1) # increment infected count
        newR = np.where(I == pc.infLT + 2)
        R[newR] = 1
        I[newR] = 0
        if settings.plot:
            if t % settings.pltFreq == 0:
                pltSim(S=S, I=I, R=R, t=t+1, anim=settings.anim, show=settings.show)

    # clean metrics
    metrics.endT = t
    metrics.numS = metrics.numS[:t]
    metrics.numI = metrics.numI[:t]
    metrics.numR = metrics.numR[:t]
    metrics.maxD = metrics.maxD[:t]
    metrics.R0 = metrics.numI[1:]/metrics.numI[:-1]
    if pc.rho == 0:
        metrics.mortality_ratio = 0
    else:
        metrics.mortality_ratio = (metrics.numR[-1] + metrics.numI[-1]) / (pc.rho * pc.L**2)

    if settings.plot:
        pltSim(S=S, I=I, R=R, t=t, anim=settings.anim, show=settings.show)
    return [S,I,R], metrics
