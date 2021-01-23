"""
Methods used in the model and ensemble-averaging.
"""
import math
import numpy as np

ijDistance = lambda i, j : np.sqrt((i[0] - j[0])**2 + (i[1] - j[1])**2)  # where i and j are tuples
modelSelector = lambda m, exponent : exponent if m == 'exp' else 0.5*exponent**2


def setFields(L:int, rho:float, epiC:int, r:int, init_n_infected:int, epiType='centralised') -> tuple:
    """
    Initialise the domain in terms of three fields, S,I and R.
    """
    S = np.zeros(shape=[L, L])  # Init domain with density rho
    S_tree_number = rho * L ** 2
    tree = 0
    while tree < S_tree_number:  # seed exact number in random locations
        rand_row = np.random.randint(0, L - 1)
        rand_col = np.random.randint(0, L - 1)
        assert rand_row or rand_col < L
        if not S[rand_row, rand_col]:
            S[rand_row, rand_col] = 1
            tree += 1

    I = np.zeros_like(S)  # Init Infected field
    if epiType == 'centralised':
        if r == 0:
            I[epiC, epiC] = 1
        else:
            I[epiC - r:epiC + r, epiC - r:epiC + r] = 2
        S[np.where(I)] = 0

    elif epiType == 'distributed':
        randRow = np.array(np.random.randint(1, I.shape[0], size=init_n_infected))
        randCol = np.array(np.random.randint(1, I.shape[0], size=init_n_infected))
        randInf = tuple([randRow, randCol])
        I[randInf] = 2
        S[randInf] = 0

    R = np.zeros_like(I)  # Init removed field
    return S,I,R


def set_R0trace(I_ind: np.array, R0_hist:dict) -> None:
    for i in range(len(I_ind[0])):
        site = str(I_ind[0][i]) + str(I_ind[1][i])
        R0_hist[site] = [0, 0]


def update_R0trace(R0_hist:dict, new_trace:list, site: list) -> int:
    """
    Update the record of secondary infections, ie { str(site) : [R0, generation] }
    :param R0_trace:
    :param new_trace:
    :param site:
    :return:
    """
    infected_site_key = str(site[0]) + str(site[1])
    R0_hist[infected_site_key][0] += len(new_trace[0])  # update the source infected R0 count
    gen = R0_hist[infected_site_key][1] + 1
    for i in range(len(new_trace[0])):
        # initialise new, n^th order, infection into the record
        R0_hist[str(new_trace[0][i]) + str(new_trace[1][i])] = [0, gen]
    return gen


def ith_new_infections(infected_site: list, S_ind:np.array, alpha:float, model:str, ell:float, beta:float) -> np.array:
    """
    Get the new infections due to the ith infected tree
    :param infected_site:
    :param S_ind: Susceptible tree indices
    :param alpha: lattice constant
    :param model: type of dispersal
    :param ell: dispersal constant
    :param beta: infectivity parameter
    :return:
    """
    distance = ijDistance(i=infected_site, j=S_ind) * alpha  # (m) get distance of infected to all S-trees
    exponent = modelSelector(m=model, exponent=(distance / ell))  # d/ell (dimensionless)
    prS_I = np.exp(-exponent) * beta  # Gaussian or exponential dispersal
    ijInfected = np.array(np.random.uniform(low=0, high=1, size=len(prS_I)) < prS_I).astype(int)
    return np.where(ijInfected)

def get_new_I(S_ind: np.array, I_ind:np.array, beta:float, alpha:float, ell:float,
              model:str, R0_histories:dict) -> tuple:
    """
    Return a list of indices of all the newly infected trees, along with the max infectious order
    """
    R0_count = 0
    S_t0 = len(S_ind[0])
    newI_ind = [[],[]]
    max_gen = 0
    for i in range(len(I_ind[0])):
        # for each infected site, find secondary infections
        infected_site = [I_ind[0][i], I_ind[1][i]]
        new_I = ith_new_infections(infected_site, S_ind, alpha, model, ell, beta)
        newI_ind[0].extend(S_ind[0][new_I])  # extend the newly infected list
        newI_ind[1].extend(S_ind[1][new_I])
        gen = update_R0trace(R0_histories, new_trace=[S_ind[0][new_I], S_ind[1][new_I]], site=infected_site)
        if gen > max_gen:
            max_gen = gen
        S_ind = tuple([np.delete(S_ind[0], new_I), np.delete(S_ind[1], new_I)])
        R0_count += len(new_I[0])

    assert R0_count == S_t0 - len(S_ind[0])
    return tuple(newI_ind), max_gen


def timerPrint(seconds: float) -> str:
    """
    Print the number of hours, minutes and seconds it took to run a simulation, or ensemble.
    """
    dc = round(seconds - math.floor(seconds), 2)
    seconds = math.floor(seconds)
    hrs = math.floor(seconds / 3600)
    mns = math.floor(seconds / 60)
    secs = seconds % 60
    if seconds < 60:
        return f'{round(secs + dc, 2)} (s)'
    elif 60 <= seconds < 3600:
        return f'{mns} (mins): {secs} (s)'
    elif seconds >= 3600:
        return f'{hrs} (Hrs): {mns%60} (mins): {secs} (s)'


def R0_generation_mean(R0_trace: dict) -> list:
    """
    From the infectious history of all infected trees, calculate the generational mean.
    """
    R0_count = [0 for i in range(1000)]
    num_trees_in_gen = [0 for i in range(1000)]
    max_gen_in_sim = 0
    for site in R0_trace:
        inf_hist = R0_trace[site] # inf_hist[0]: the number of secondary infections, inf_hist[1]: the generation
        R0_count[inf_hist[1]] += inf_hist[0]
        num_trees_in_gen[inf_hist[1]] += 1  # cumulatively add the number of trees in each infectious generation
        if inf_hist[1] > max_gen_in_sim:
            max_gen_in_sim = inf_hist[1]  # update the highest generation

    R0_count = R0_count[:max_gen_in_sim]
    num_trees_in_gen = num_trees_in_gen[:max_gen_in_sim]
    return [R0_c/num_trees for R0_c, num_trees in zip(R0_count, num_trees_in_gen)]  # Return mean infections / gen


def avg_multi_dim(arrays: np.array) -> np.array:
    """
    Average arrays with variable entries.
    :param arrays: variable number of array-like structs
    """
    max_= max([len(arr) for arr in arrays])
    avg_arr = np.zeros(max_)
    for arr in arrays:
        avg_arr[0:len(arr)] += arr
    return avg_arr/len(arrays)

