"""
Define classes for model parameters, settings and metrics.
"""
import os
import numpy as np

PATH_TO_DATA_STORE = f'{os.getcwd()}/temp_dat_store/'


class ModelParamSet:  # Set model parameter and simulation fields
    META_DATA = None
    alpha = 5  # (m) lattice scale parameter
    L = 1000  # L x L = domain dim : modelled size = alpha^2 * L^2 (m^2)
    infected_lt = 50  # (steps) infectious life-time
    tend = 100  # (steps) final end time
    betas = [0.001]  # (step^-1) infectivity parameter
    rhos = [0.01]   # tree density
    ell = 1000  # (m) dispersal
    epi_center = int(L/2)
    init_n_infected = 1  # control randomly distributed epicenters
    r = 0  # control centralised epicenter radius
    model = ['exponential', 'gaussian', 'power_law'][1]


class Metrics:
    endT = 0
    max_gen = None

    save_end_time = True
    save_percolation = True
    save_mortality_ratio = True

    track_R0_history = True
    track_time_series = True

    # self._R0 = np.zeros(3000)
    # self._R0_histories = {}
    # self._maxD = np.zeros(3000)
    # self._numI = np.zeros(3000)
    # self._numS = np.zeros(3000)
    # self._numR = np.zeros(3000)


class Settings:  # Simulation setup
    plot = False
    show = False
    anim = False
    pltFreq = 2
    verbose = 3  # verbosity
    ext = '.png'
    percolation_bcd = False  # Terminate upon infection reaching the boundary
    gen_limit = 1  # if no infected trees of order `gen-limit' are  left, end the simulation.


ParamsAndSetup = {'params': ModelParamSet, 'setup': Settings, 'metrics': Metrics}




