import numpy as np
from typing import Union
from tree_epi_dispersal.plot_methods import plt_sim_frame
from parameters_and_settings import ModelParamSet, Settings, Metrics
from tree_epi_dispersal.model_dynamics_helpers import set_R0trace, ijDistance, setFields, ith_new_infections

printStep = lambda t, freq : print('\t\t Time : {} (days)'.format(t)) if t % freq == 0 else None


def get_new_I(S_ind: np.array, I_ind: np.array, beta: float, ell: float, R0_histories: dict) -> tuple:
    """
    Return a list of indices of all the newly infected trees, along with the max infectious order
    """
    alpha = ModelParamSet.alpha
    model = ModelParamSet.model
    R0_count = 0
    num_S = len(S_ind[0])
    newI_ind = [[], []]
    max_generation = Settings.max_generation_bcd if R0_histories else None
    max_gen_exceeded = False if R0_histories else None

    # the maximum generation of infected trees considered
    for i in range(len(I_ind[0])):
        # for each infected site, find secondary infections
        infected_site = [I_ind[0][i], I_ind[1][i]]
        new_I = ith_new_infections(infected_site, S_ind, alpha, model, ell, beta)
        newI_ind[0].extend(S_ind[0][new_I])  # extend the newly infected list
        newI_ind[1].extend(S_ind[1][new_I])

        if R0_histories:
            print('tracking R0...')
            #todo sort me out.....
            assert 0
            gen = update_R0trace(R0_histories, new_trace=[S_ind[0][new_I], S_ind[1][new_I]], site=infected_site)
            if gen_limit is not None and gen <= max_generation:
                max_gen_exceeded = False  # if a single tree, of less than or equal to, order gen exists continue simulation
            continue

        S_ind = tuple([np.delete(S_ind[0], new_I), np.delete(S_ind[1], new_I)])
        R0_count += len(new_I[0])

    assert R0_count == num_S - len(S_ind[0])
    return tuple(newI_ind), max_gen_exceeded


def run_simulation(rho: float, beta: float, ell: Union[float, tuple]) -> dict:
    """

    :param rho: tree density
    :param beta: pathogen infectiviy
    :param ell: pathogen dispersal parameter(s)
    :return: sim_result, a dictionary of required fields.
    """

    S, I, R = setFields(rho)
    R0_history = set_R0trace(I, {}) if Metrics.track_R0_history else None

    t = None
    percolation_event = True
    all_infected_trees_died = False

    S_ts = np.zeros(ModelParamSet.tend) if Metrics.track_time_series else None
    I_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None
    R_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None
    max_d_ts = np.zeros_like(S_ts) if Metrics.track_time_series else None
    epi_c = ModelParamSet.epi_center if Metrics.track_time_series else None

    for t in range(ModelParamSet.tend):
        if Settings.verb == 2:
            printStep(t, freq=1)

        S_ = np.where(S)
        I_ = np.where(I)
        R_ = np.where(R)

        num_infected = len(I_[0])

        # BCD 1, all infected trees dies/removed
        if not num_infected:
            all_infected_trees_died = True
            print('broke all trees dead @ ', t)
            break

        if Metrics.track_time_series:
            S_ts[t] = len(S_[0])
            I_ts[t] = num_infected
            R_ts[t] = len(R_[0])
            max_d_ts[t] = ijDistance(i=[epi_c, epi_c], j=I_).max() * ModelParamSet.alpha

        if Settings.percolation_bcd and max_d_ts[t] >= (ModelParamSet.L/2 - 10) * ModelParamSet.alpha:
            if Metrics.save_percolation:
                percolation_event = True
            print('broke percolation ')
            break

        # update fields S, I, R
        newI_ind, max_gen_exceeded = get_new_I(S_, I_, beta, ell, R0_history)

        if Settings.max_generation_bcd and max_gen_exceeded:
            # if no remaining infected trees of order `gen-limit', terminate simulation
            break

        S[newI_ind] = 0
        I[newI_ind] = 1
        I = I + np.array(I >= 1).astype(int)  # increment infected count
        # Life-time dynamics : Pr I -> R = 1-exp^(-1 t/mu)
        newR = np.exp(-I/ModelParamSet.infected_lt) < np.random.uniform(0, 1, size=S.shape)
        newR = np.where(newR)
        R[newR] = 1
        I[newR] = 0

        if Settings.plot and t % Settings.plt_freq == 0:
            plt_sim_frame(S, I, R, t+1, Settings.save, Settings.show)

    # save required fields
    sim_result = {}

    if Metrics.track_R0_history:
        sim_result['R0_hist'] = R0_history

    if Metrics.track_time_series:  # clean metrics
        S_ts = S_ts[:t]
        I_ts = I_ts[:t]
        R_ts = R_ts[:t]
        max_d_ts = max_d_ts[:t]
        sim_result['time_series'] = {'S': S_ts, 'I': I_ts, 'R': R_ts, 'max_d': max_d_ts}

    if Metrics.save_end_time and all_infected_trees_died:
        sim_result['sim_end_time'] = t

    if Metrics.save_mortality_ratio:
        sim_result['mortality_ratio'] = (R_ts[-1] + I_ts[-1]) / (rho * ModelParamSet.L**2) if rho > 0 else 0

    if Settings.plot:
        plt_sim_frame(S, I, R, t, Settings.save, Settings.show)

    if Metrics.save_percolation:
        sim_result['percolation'] = percolation_event

    return sim_result
