import numpy as np

ijDistance = lambda i, j : np.sqrt((i[0] - j[0])**2 + (i[1] - j[1])**2)
modelSelector = lambda m, exponent : exponent if m == 'exp' else 0.5*exponent**2
printStep = lambda t, freq : print('\t Time : {} (days)'.format(t)) if t % freq == 0 else None

def prDispersal(S_ind, I_ind, alpha, ell, L, model) -> 'Pr(S_x-> I_x| I_y)':
    """
    :return: [L, L] float-type array of dispersal probabilities
    """
    prS_S = np.ones_like(S_ind[0])
    for i in range(len(I_ind[0])):
        # for each infected site, find Pr of infecting S neighbours
        infected_site = [I_ind[0][i], I_ind[1][i]]
        distance = ijDistance(i=infected_site, j=S_ind) * alpha  # (m)
        exponent = modelSelector(m=model, exponent=(distance/ell))  # d/ell (dimensionless)
        prS_I = np.exp(-exponent)  # Gaussian or exponential dispersal
        prS_S = prS_S * ( 1 - prS_I)

    prS_I = 1 - prS_S
    prDisp = np.zeros(shape=(L, L)).astype(float)
    prDisp[S_ind] = prS_I
    return prDisp

def findR0(pc, numR, settings) -> 'R0_max':
    """
    Simulate number of secondary infections due to one infectious case
    :return:
    """
    if settings.plot:
        from plots.plotLib import pltSim
    for t in range(pc.infLT):
        S_ = np.where(pc.S)
        I_ = np.where(pc.I)
        R_ = np.where(pc.R)
        # update metrics
        numR[t] = len(R_[0])
        prS_I = pc.beta * prDispersal(S_, I_, pc.alpha, pc.ell, pc.L, pc.model)
        newR = pc.randGen(dim=pc.L, comp=prS_I)
        pc.S[np.where(newR)] = 0
        pc.R[np.where(newR)] = 1
        # plot simulation fields
        if settings.plot:
            if t % settings.pltFreq == 0 or t==pc.infLT-1:
                pltSim(S=pc.S, I=pc.I, R=pc.R, t=t, anim=settings.anim)

    # metrics.endT = t
    # metrics.numR = metrics.numR[:t]
    return numR

return_ = {'pc': 'with updated output fields', 'metrics': 'holding time-series data'}
def runSim(pc, metrics, settings) -> return_:
    """
    Run dynamic simulation of pathogen spread
    :return:
    """
    if settings.plot:
        from plots.plotLib import pltSim
    pc.setField
    for t in range(pc.tend):
        if settings.verbose2:
            printStep(t, freq=10)
        S_ = np.where(pc.S)
        I_ = np.where(pc.I)
        R_ = np.where(pc.R)
        # update metrics
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
        prS_I = pc.beta * prDispersal(S_, I_, pc.alpha, pc.ell, pc.L, pc.model)
        newI = pc.randGen(L=pc.L, thresh=prS_I)
        pc.S[np.where(newI)] = 0
        pc.I = pc.I + (pc.I > 1) + newI * 2  # increment infected count
        newR = np.where(pc.I == pc.infLT + 1)
        pc.R[newR] = 1
        pc.I[newR] = 0
        # plot simulation fields
        if settings.plot:
            if t % settings.pltFreq == 0:
                pltSim(S=pc.S, I=pc.I, R=pc.R, t=t, anim=settings.anim, show=settings.show)

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


if __name__ == '__main__':
    # test model
    import matplotlib.pyplot as plt
    ell = 25
    alpha = 5
    S = [np.linspace(0, 20, 100), np.zeros(100)]
    distance = ijDistance(j=S, i=[0, 0]) * alpha
    plt.plot(distance, np.exp(-distance/ell), label='exp dist')
    plt.plot([ell, ell], [0, 1/np.e])
    plt.plot(distance, np.exp(-distance**2/ell**2), label='normal dist')
    plt.legend()
    plt.show()