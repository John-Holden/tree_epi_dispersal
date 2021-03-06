import numpy as np
from tree_epi_dispersal.model_dynamics_helpers import set_R0trace, ijDistance, get_new_I, setFields
from parameters_and_settings import ModelParamSet, Settings, Metrics

printStep = lambda t, freq : print('\t\t Time : {} (days)'.format(t)) if t % freq == 0 else None


def runSim(rho: float, beta: float) -> '[S,I,R], metrics':
    """
    Run dynamic simulation of pathogen spread,
    :return: SIR fields and updated metrics
    """

    S, I, R = setFields(rho)
    R0_history = set_R0trace(I, {}) if Metrics.track_R0_history else None
    S_ts = np.zeros(ModelParamSet.tend) if Metrics.track_time_series else None
    I_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None
    R_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None
    epi_c = ModelParamSet.epi_center if Metrics.track_time_series else None
    max_d_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None

    if Settings.plot:
        from tree_epi_dispersal.plot_methods import pltSim

    for t in range(ModelParamSet.tend):
        if Settings.verbose == 2:
            printStep(t, freq=1)
        S_ = np.where(S)
        I_ = np.where(I)
        R_ = np.where(R)
        # update metrics
        num_infected = len(I_[0])

        if Metrics.track_time_series:
            S_ts = len(S_[0])
            I_ts = num_infected
            R_ts = len(R_[0])
            max_d_ts[t] = ijDistance(i=[epi_c, epi_c], j=I_).max() * ModelParamSet.alpha

        if not num_infected:  # BCD 1, all infected trees dies/removed
            Metrics.all_infected_trees_died = True
            Metrics.simulation_end_time = t if Metrics.save_end_time else None
            break

        if Settings.percolation_bcd and max_d_ts[t] >= (ModelParamSet.L/2 - 10) * ModelParamSet.alpha:
            if Metrics.save_percolation:
                Metrics.time_percolated = t
                Metrics.percolation_event = True
            break

        # update fields S, I, R
        newI_ind, max_gen_exceeded = get_new_I(S_, I_, pc.beta, pc.alpha, pc.ell, pc.model,
                                                   metrics._R0_histories, Settings.gen_limit)

        if Settings.gen_limit is not None and max_gen_exceeded:
            # if no remaining infected trees of order `gen-limit', terminate simulation
            break

        S[newI_ind] = 0
        I[newI_ind] = 1
        I = I + (I>=1) # increment infected count
        newR = np.exp(-I/pc.infLT) < np.random.uniform(0, 1, size=S.shape)  # Pr I -> R = 1-exp^(-1 t/mu)
        newR = np.where(newR)
        R[newR] = 1
        I[newR] = 0
        if Settings.plot:
            if t % Settings.pltFreq == 0:
                pltSim(S=S, I=I, R=R, t=t+1, anim=Settings.anim, show=Settings.show)

    # clean metrics
    metrics.endT = t
    metrics._numS = metrics._numS[:t]
    metrics._numI = metrics._numI[:t]
    metrics._numR = metrics._numR[:t]
    metrics._maxD = metrics._maxD[:t]
    metrics.R0 = metrics._numI[1:]/metrics._numI[:-1]
    if pc.rho == 0:
        metrics.mortality_ratio = 0
    else:
        metrics.mortality_ratio = (metrics._numR[-1] + metrics._numI[-1]) / (pc.rho * pc.L**2)

    if Settings.plot:
        pltSim(S=S, I=I, R=R, t=t, anim=Settings.anim, show=Settings.show)

    return [S,I,R], metrics
