"""
Define classes for model parameters, settings and metrics.
"""
import os
import numpy as np

PATH_TO_DATA_STORE = f'{os.getcwd()}/temp_dat_store/'


class ModelParamSet:  # Set model parameter and simulation fields
    META_DATA = 'Exponentially distributed life-time dynamics : True'
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


class Metrics:   # Define which metrics are recorded over the simulation
    max_gen = None
    save_end_time = True
    save_percolation = False
    save_mortality_ratio = True
    track_R0_history = False
    track_time_series = True


class Settings:  # Simulation setup
    plot = False
    show = False
    save = False
    plt_freq = 2
    verb = 2  # verbosity
    ext = '.png'
    percolation_bcd = False  # Terminate upon infection reaching the boundary
    max_generation_bcd = False  # if no infected trees of order @ max_gen are  left, end the simulation.


ParamsAndSetup = {'params': ModelParamSet, 'setup': Settings, 'metrics': Metrics}




