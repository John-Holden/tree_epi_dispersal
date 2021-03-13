"""
Methods used in the model and ensemble-averaging.
"""

import numpy as np
from typing import Callable
from parameters_and_settings import ModelParamSet


def ij_distance(i: tuple, j: tuple) -> np.ndarray:
    """
    Find the distance between site i : (x, y) and coordinates j : (x1,...xN), (y1,...,yN)
    """
    return np.sqrt((i[0] - j[0]) ** 2 + (i[1] - j[1]) ** 2)


def ij_arr_distance(xy: tuple, arr: np.ndarray) -> np.ndarray:
    """
    Compute the distance of points p in an array to the origin xy
    :param xy: position from which distance matrix will be calculated
    :param arr: any array
    :return:
    """
    xarr, yarr = np.meshgrid(np.arange(0, arr.shape[0], 1),
                             np.arange(0, arr.shape[1], 1))
    return np.sqrt((xarr - xy[0])**2 + (yarr - xy[1])**2)


def model_selector() -> Callable:
    """
    Define the current dispersal model
    :return: dispersal function
    """
    model = None

    ModelParamSet.assert_config()

    if ModelParamSet.model == 'gaussian':
        def model(dist, beta, ell):
            return beta * np.exp(-0.5*(dist/ell) ** 2)
    elif ModelParamSet.model == 'exponential':
        def model(dist, beta, ell):
            return beta * np.exp(-(dist/ell))
    elif ModelParamSet.model == 'power_law':
        def model(dist, beta, ell):
            a, b = ell
            return beta * (1 + dist/a)**(-b)

    return model


def set_SIR(S: np.ndarray, epicenter_init_cond: str):
    """
    Set the fields used for the SIR model
    :param S: susceptible tree's
    :param epicenter_init_cond: initial conditions of epicenters
    :return:
    """
    I = np.zeros_like(S)  # Infected field
    r = ModelParamSet.r
    if epicenter_init_cond == 'centralised':
        epi_c = ModelParamSet.epi_center
        if r == 0:
            I[epi_c, epi_c] = 1
        else:
            I[epi_c - r:epi_c + r, epi_c - r:epi_c + r] = 2
        S[np.where(I)] = 0

    elif epicenter_init_cond == 'distributed':
        num_init_infected = ModelParamSet.init_n_infected
        randRow = np.array(np.random.randint(1, I.shape[0], size=num_init_infected))
        randCol = np.array(np.random.randint(1, I.shape[0], size=num_init_infected))
        randInf = tuple([randRow, randCol])
        I[randInf] = 1
        S[randInf] = 0

    R = np.zeros_like(I)  # Init removed field
    return S, I, R


def set_ADB(S_tr: np.ndarray, max_epi: int = 10) -> 'tuple of np.array-fields':
    """
    Set the fields used for the ash dieback model.
    This includes a source-infected tree surrounded by a small number of infectious fruiting bodies.
    A random number of fruiting body-sources, between 1-N, is sampled from a uniform distribution and
    placed in the neighbourhood of the infectious tree.
    :param S_tr: susceptible tree's
    :return:
    """
    I_fb = np.zeros_like(S_tr)  # Infectious fruiting bodies

    epi_c = ModelParamSet.epi_center
    dist_pr = np.exp(-ij_arr_distance((epi_c, epi_c), S_tr) / ModelParamSet.alpha*2)
    set_epi_c, epi_n = 0, np.random.randint(1, max_epi+1)

    potential_epi_c = dist_pr > np.random.uniform(low=0, high=1, size=S_tr.shape)

    potential_epi_c = np.where(potential_epi_c)
    num_potential = len(potential_epi_c[0])
    assert num_potential > epi_n, f'Error, potential epicenter number {len(potential_epi_c[0])} < {epi_n}'

    while set_epi_c < epi_n:
        row, col = np.random.randint(0, num_potential), np.random.randint(0, num_potential)
        row, col = potential_epi_c[0][row], potential_epi_c[1][col]

        if I_fb[row, col]:
            continue

        I_fb[row, col] = 1
        set_epi_c += 1

    S_tr[epi_c, epi_c] = 0
    S_tr = np.where(S_tr)
    I_fb = np.where(I_fb)

    # define life-times of fruiting bodies
    I_fb_lf = np.random.normal(loc=ModelParamSet.fb_lt[0], scale=ModelParamSet.fb_lt[1], size=len(I_fb[0])).astype(int)

    I_fb = [I_fb[0], I_fb[1], I_fb_lf]
    I_tr = ([epi_c], [epi_c])  # Infectious trees
    E_tr = ([], [])  # Exposed/latently infected hosts
    R_fb = [[], [], []]  # Removed fruiting body sources

    return S_tr, I_tr, E_tr, I_fb, R_fb


def setFields(rho: float, model: str, epi_IC='centralised') -> tuple:
    """
    Initialise the domain-fields
    """
    S = np.zeros(shape=[ModelParamSet.L, ModelParamSet.L])  # init susceptible domain with density rho
    S_tree_number = rho * ModelParamSet.L ** 2
    tree = 0

    while tree < S_tree_number:  # seed exact number in random locations
        rand_row = np.random.randint(0, ModelParamSet.L)
        rand_col = np.random.randint(0, ModelParamSet.L)

        if not S[rand_row, rand_col]:
            S[rand_row, rand_col] = 1
            tree += 1

    if model == 'SIR':
        return set_SIR(S, epi_IC)

    elif model == 'ADB':
        return set_ADB(S)

    else:
        raise ValueError('Error, wrong input')


def set_R0trace(I: np.array, R0_hist: dict, test_mode: bool, adb_mode: bool = False) -> dict:
    """
    Initialise the R0-history dictionary -- record the generation of infect along with infection statistics
        - site_ij : [generation_infected, infectious_count, [distance_of_infected_trees, `optional`]]
    """
    I_ind = np.where(I) if not adb_mode else I
    for i in range(len(I_ind[0])):
        site = str(I_ind[0][i]) + str(I_ind[1][i])
        R0_hist[site] = [0, 0, []] if test_mode else [0, 0]

    return R0_hist


def update_R0trace(R0_hist: dict, new_infected: tuple, site: tuple, test_mode: bool,
                   update_secondaries: bool = True) -> int:
    """
    Update the record of secondary infections, i.e:
        - site_i_j : [R0, generation]
        - site_i_j is a string-representation of coordinates
        - R0 is the number of current total number of infections
        - generation is the n^th order generation site_i_j became infected

    :param R0_hist: store infection counts due to each infected, along with generation became infected
    :param new_infected: indices of newly infected trees
    :param site: coordinates (i, j) of source-infection
    :param update_secondaries: if true, update initialise new secondary infections
    :return:
    """

    source_infection = str(site[0]) + str(site[1])
    R0_hist[source_infection][0] += len(new_infected[0])  # update the source infected R0 count
    if test_mode:
        # append extra information about the dispersal
        R0_hist[source_infection][2].extend(ij_distance(site, new_infected) * ModelParamSet.alpha)

    generation_of_new_infection = R0_hist[source_infection][1] + 1
    if update_secondaries:
        for i in range(len(new_infected[0])):
            # initialise new, n^th order, infection into the record
            R0_hist[str(new_infected[0][i]) + str(new_infected[1][i])] = [0, generation_of_new_infection]

    return generation_of_new_infection


def ith_new_infections(infected_site: tuple, S_ind: np.array, ell: float, beta: float,
                       dispersal_model: Callable) -> np.array:
    """
    Get the new infections due to the ith infected tree
    :param infected_site: coordinates (i, j) of infected site
    :param S_ind: Susceptible tree indices
    :param dispersal_model: type of dispersal
    :param ell: dispersal constant(s)
    :param beta: infectivity parameter
    :return: indices of newly infected trees
    """
    alpha = ModelParamSet.alpha
    # (m) get distance of infected to all S-trees
    distance = ij_distance(i=infected_site, j=S_ind) * alpha
    # probability of susceptible trees S transitioning to infected
    pr_S_I = dispersal_model(distance, beta, ell)
    ijInfected = np.array(np.random.uniform(low=0, high=1, size=len(pr_S_I)) < pr_S_I).astype(int)
    return np.where(ijInfected)


def R0_generation_mean(R0_trace: dict) -> list:
    """
    From the infectious history of all infected trees, calculate the mean for each generation i.
    """
    R0_count = np.zeros(1000)
    num_trees_in_gen = np.zeros(1000)
    max_gen_in_sim = 0
    for site in R0_trace:
        # inf_hist[0]: the number of secondary infections, inf_hist[1]: the generation
        inf_hist = R0_trace[site]
        R0_count[inf_hist[1]] += inf_hist[0]
        num_trees_in_gen[inf_hist[1]] += 1  # cumulatively add the number of trees in each infectious generation

        max_gen_in_sim = inf_hist[1] if inf_hist[1] > max_gen_in_sim else max_gen_in_sim  # update the max generation

    R0_count = R0_count[:max_gen_in_sim]
    num_trees_in_gen = num_trees_in_gen[:max_gen_in_sim]
    return [R0_c / num_trees for R0_c, num_trees in zip(R0_count, num_trees_in_gen)]  # Return mean infections / gen


def avg_multi_dim(arrays: np.array) -> np.array:
    """
    Average arrays with variable entries.
    :param arrays: variable number of array-like structs
    """
    max_ = max([len(arr) for arr in arrays])
    avg_arr = np.zeros(max_)
    for arr in arrays:
        avg_arr[0:len(arr)] += arr
    return avg_arr / len(arrays)
