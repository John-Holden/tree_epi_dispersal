import numpy as np
import datetime
from typing import Union, Callable
from tree_epi_dispersal.plot_methods import plt_sim_frame, plt_adb_frame
from parameters_and_settings import ModelParamSet, Settings, Metrics
from tree_epi_dispersal.model_dynamics_helpers import set_R0trace, ij_distance, setFields, ith_new_infections, \
    update_R0trace, model_selector

printStep = lambda t, freq: print('\t\t Time : {} (days)'.format(t)) if t % freq == 0 else None


def get_new_I(S_ind: np.array, I_ind: np.array, beta: float, ell: float, R0_histories: dict,
              dispersal_model: Callable, test_mode: bool, update_secondaries: bool = True) -> tuple:
    """
    Return a list of indices of all the newly infected trees, along with the max infectious order
    """
    R0_count = 0
    num_S = len(S_ind[0])
    newI_ind = [[], []]
    max_gen_exceeded = True if Settings.max_generation_bcd else None
    # the maximum generation of infected trees considered
    for i in range(len(I_ind[0])):
        # for each infected site, find secondary infections
        infected_site = (I_ind[0][i], I_ind[1][i])
        new_I = ith_new_infections(infected_site, S_ind, ell, beta, dispersal_model)
        newI_ind[0].extend(S_ind[0][new_I])  # extend the newly infected list
        newI_ind[1].extend(S_ind[1][new_I])
        if R0_histories:
            S_trans_I = (S_ind[0][new_I], S_ind[1][new_I])  # sites S that will transition to the infected-state
            generation_of_infected_site = update_R0trace(R0_histories, S_trans_I, infected_site,
                                                         test_mode, update_secondaries)
            # if no trees of `max_Gen` are left, terminate simulation
            if Settings.max_generation_bcd and generation_of_infected_site <= Settings.max_generation_bcd:
                max_gen_exceeded = False

        S_ind = tuple([np.delete(S_ind[0], new_I), np.delete(S_ind[1], new_I)])  # take away newly infected from S
        R0_count += len(new_I[0])

    assert R0_count == num_S - len(S_ind[0])

    return S_ind, tuple(newI_ind), max_gen_exceeded


def run_ADB(rho: float, beta: float, ell: Union[int, float, tuple]):
    """
    Simulate a configured ash dieback dispersal_model
    :param rho: tree density
    :param beta: a compound parameter representing pathogen infectiousness
    :param ell: pathogen dispersal_type parameter(s)
    :return: sim_result, a dictionary of required fields.
    """

    S_tr, I_tr, E_tr, I_fb, R_fb = setFields(rho, model='ADB', epi_IC='distributed')
    # track the number of infections due to the fruiting bodies
    R0_history = set_R0trace(I_fb, {}, test_mode=False, adb_mode=True)
    dispersal_model = model_selector()  # get function for the current kernel for the configuration
    break_condition = None
    start_date = datetime.datetime(2020, 6, 1)

    for t in range(ModelParamSet.tend):
        current_date = start_date + datetime.timedelta(days=t)
        if Settings.verb == 2:
            print('t : ', current_date.strftime("%b %d"))

        S_tr, new_E_tr, max_gen_exceeded = get_new_I(S_tr, I_fb, beta, ell, R0_history, dispersal_model,
                                                     test_mode=False, update_secondaries=False)

        # S_tree -> E_tree
        E_tr[0].extend(new_E_tr[0]), E_tr[1].extend(new_E_tr[1])

        # I_fruiting_body -> R_fruiting_body
        to_remove_fb = [index for index, to_rem in enumerate(I_fb[2] < t) if to_rem]
        R_fb[0].extend([I_fb[0][ind] for ind in to_remove_fb])
        R_fb[1].extend([I_fb[1][ind] for ind in to_remove_fb])
        R_fb[2].extend([t - 1 for ind in to_remove_fb])

        I_fb[0] = np.delete(I_fb[0], to_remove_fb)
        I_fb[1] = np.delete(I_fb[1], to_remove_fb)
        I_fb[2] = np.delete(I_fb[2], to_remove_fb)

        if not len(I_fb[0]):
            break_condition = 'fruiting body extinction'
            break

        if Settings.plot and t % Settings.plot_freq == 0:
            plt_adb_frame(S_tr, E_tr, I_fb, R_fb, current_date.strftime("%b %d"))

    sim_result = {"termination": break_condition,
                  'R0_hist': R0_history}

    return sim_result


def run_SIR(rho: float, beta: float, ell: Union[int, float, tuple],
            test_mode: bool = False) -> dict:
    """
    Simulate a generic SIR-dispersal_type dispersal_model
    :param rho: tree density
    :param beta: a compound parameter representing pathogen infectiousness
    :param ell: pathogen dispersal_type parameter(s)
    :param test_mode : if true, limit the dynamics so only the initially infected infect neighbours
    :return: sim_result, a dictionary of required fields.
    """

    S, I, R = setFields(rho, model='SIR')
    R0_history = set_R0trace(I, {}, test_mode) if Metrics.save_R0_history else None
    t = None
    percolation_event = True
    all_infected_trees_died = False

    S_ts, I_ts, R_ts, max_d_ts, epi_c = None, None, None, None, None
    if Metrics.save_time_series or Metrics.save_mortality_ratio:
        S_ts = np.zeros(ModelParamSet.tend)
        I_ts = np.zeros(ModelParamSet.tend)
        R_ts = np.zeros(ModelParamSet.tend)
        max_d_ts = np.zeros_like(S_ts)
        epi_c = ModelParamSet.epi_center

    model = model_selector()  # get function for the current kernel for the configuration
    break_condition = None

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
            break_condition = 'bcd1 : all trees dead'
            break

        if Metrics.save_time_series or Metrics.save_mortality_ratio:
            S_ts[t] = len(S_[0])
            I_ts[t] = num_infected
            R_ts[t] = len(R_[0])
            max_d_ts[t] = ij_distance(i=(epi_c, epi_c), j=I_).max() * ModelParamSet.alpha

        if Settings.percolation_bcd and max_d_ts[t] >= (ModelParamSet.L/2 - 10) * ModelParamSet.alpha:
            if Metrics.save_percolation:
                percolation_event = True
            break_condition = 'bcd2: percolation '
            break

        # update fields S, I, R
        _, newI_ind, max_gen_exceeded = get_new_I(S_, I_, beta, ell, R0_history, model, test_mode)
        newI_ind = ([], []) if test_mode else newI_ind  # if test mode prevent secondary infections from infecting

        if Settings.max_generation_bcd and max_gen_exceeded:
            # if no remaining infected trees of order `gen-limit', terminate simulation
            break_condition = f'bcd3: no more trees of gen {Settings.max_generation_bcd} left'
            break

        S[newI_ind] = 0
        I[newI_ind] = 1
        I = I + np.array(I >= 1).astype(int)  # increment infected count
        # Life-time dynamics : Pr I -> R = 1-exp^(-1 t/mu)
        newR = np.exp(-I/ModelParamSet.infected_lt) < np.random.uniform(0, 1, size=S.shape)
        newR = np.where(newR)
        R[newR] = 1
        I[newR] = 0

        if Settings.plot and t % Settings.plot_freq == 0:
            plt_sim_frame(S, I, R, t+1, Settings.save, Settings.show)

    break_condition = 'bcd0: complete simulation time elapsed' if break_condition is None else break_condition
    # save required fields
    sim_result = {"termination": break_condition}

    if Metrics.save_R0_history:
        sim_result['R0_hist'] = R0_history

    if Metrics.save_time_series:  # clean metrics
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
        plt_sim_frame(S, I, R, t, Settings.save, Settings.show, msg=' : Out')

    if Metrics.save_percolation:
        sim_result['percolation'] = percolation_event

    return sim_result
