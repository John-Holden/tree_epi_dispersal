import numpy as np
import sys

ijDistance = lambda i, j : np.sqrt((i[0] - j[0])**2 + (i[1] - j[1])**2)
modelSelector = lambda m, exponent : exponent if m == 'exp' else 0.5*exponent**2
printStep = lambda t, freq : print('\t Time : {} (days)'.format(t)) if t % freq == 0 else None


def set_R0trace(I_ind, R0_trace):
    for i in range(len(I_ind[0])):
        site = str(I_ind[0][i]) + str(I_ind[1][i])
        R0_trace[site] = [0, 0]
    return

def update_R0trace(R0_trace, new_trace, site) -> '{ [ site ] : [R0, generation]}':
    R0_trace[str(site[0]) + str(site[1])][0] += len(new_trace[0])
    gen = R0_trace[str(site[0]) + str(site[1])][1] + 1
    for i in range(len(new_trace[0])):
        # initialise nth generation infections
        R0_trace[str(new_trace[0][i]) + str(new_trace[1][i])] = [0, gen]
    return

def ith_new_infections(infected_site, S_ind, alpha, model, ell, beta) -> 'new infections due to i':
    distance = ijDistance(i=infected_site, j=S_ind) * alpha  # (m)
    exponent = modelSelector(m=model, exponent=(distance / ell))  # d/ell (dimensionless)
    prS_I = np.exp(-exponent) * beta  # Gaussian or exponential dispersal
    ijInfected = np.array(np.random.uniform(low=0, high=1, size=len(prS_I)) < prS_I).astype(int)
    return np.where(ijInfected)

def get_new_I(S_ind, I_ind, beta, alpha, ell, model, R0_trace) -> 'Indices of all newly infected':
    R0_count = 0
    S_t0 = len(S_ind[0])
    newI_ind = [[],[]]
    for i in range(len(I_ind[0])):
        # for each infected site, find secondary infections
        infected_site = [I_ind[0][i], I_ind[1][i]]
        new_I = ith_new_infections(infected_site, S_ind, alpha, model, ell, beta)
        newI_ind[0].extend(S_ind[0][new_I])
        newI_ind[1].extend(S_ind[1][new_I])
        update_R0trace(R0_trace, new_trace=[S_ind[0][new_I], S_ind[1][new_I]], site=infected_site)
        S_ind = tuple([np.delete(S_ind[0], new_I), np.delete(S_ind[1], new_I)])
        R0_count += len(new_I[0])

    assert R0_count == S_t0 - len(S_ind[0])
    return tuple(newI_ind)

return_ = {'pc': 'with updated output fields', 'metrics': 'holding time-series data'}
def runSim(pc, metrics, settings) -> return_:
    """
    Run dynamic simulation of pathogen spread
    :return:
    """
    if settings.plot:
        from plots.plotLib import pltSim

    for t in range(pc.tend):
        if settings.verbose:
            printStep(t, freq=1)
        S_ = np.where(pc.S)
        I_ = np.where(pc.I)
        R_ = np.where(pc.R)
        # update metrics
        if t == 0:
            set_R0trace(I_, metrics.R0_trace)

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
        newI_ind = get_new_I(S_, I_, pc.beta, pc.alpha, pc.ell, pc.model, metrics.R0_trace)
        pc.S[newI_ind] = 0
        pc.I[newI_ind] = 1
        pc.I = pc.I + (pc.I>=1) # increment infected count
        newR = np.where(pc.I == pc.infLT + 2)

        pc.R[newR] = 1
        pc.I[newR] = 0
        if settings.plot:
            if t % settings.pltFreq == 0:
                pltSim(S=pc.S, I=pc.I, R=pc.R, t=t+1, anim=settings.anim, show=settings.show)

    # clean metrics
    metrics.endT = t
    metrics.numS = metrics.numS[:t]
    metrics.numI = metrics.numI[:t]
    metrics.numR = metrics.numR[:t]
    metrics.maxD = metrics.maxD[:t]
    metrics.R0 = metrics.numI[1:]/metrics.numI[:-1]
    metrics.mortality_ratio = (metrics.numR[-1] + metrics.numI[-1]) / (pc.rho * pc.L**2)
    if settings.plot:
        pltSim(S=pc.S, I=pc.I, R=pc.R, t=t, anim=settings.anim, show=settings.show)
    return pc, metrics
