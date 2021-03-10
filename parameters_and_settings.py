"""
Define classes for model parameters, settings and metrics.
"""
import os
import numpy as np

PATH_TO_DATA_STORE = f'{os.getcwd()}/temp_dat_store/'


class ModelParamSet:  # Set model parameter and simulation fields
    META_DATA = 'Exponentially distributed life-time dynamics : True'
    alpha = 5  # (m) lattice scale parameter
    L = 500  # L x L = domain dim : modelled size = alpha^2 * L^2 (m^2)
    infected_lt = 500  # (steps) infectious life-time -- exponentially-distributed mean = 1/T
    tend = 500  # (steps) final end time
    betas = [0.0005]  # (step^-1) infectivity-parameter
    rhos = [0.01]   # tree density
    # dispersal values from https://doi.org/10.1093/femsec/fiy049 are:
    #   - gaussian = 195m
    #   - power-law : a = 206, b=3.3
    # ell = (205, 3.3)  # (m) dispersal
    ell = 195
    epi_center = int(L/2)
    init_n_infected = 1  # control randomly distributed epicenters
    r = 0  # control centralised epicenter radius
    model = ['exponential', 'gaussian', 'power_law'][1]


class Metrics:   # Define which metrics are recorded over the simulation
    save_end_time = True
    save_percolation = False
    save_mortality_ratio = True
    save_R0_history = True
    save_time_series = False


class Settings:  # Simulation setup
    plot = True
    show = True
    save = False
    plot_freq = 50
    verb = 2  # verbosity
    ext = '.png'
    percolation_bcd = False  # Terminate upon infection reaching the boundary
    max_generation_bcd = False  # if no infected trees of order @ max_gen are  left, end the simulation.


ParamsAndSetup = {'params': ModelParamSet, 'setup': Settings, 'metrics': Metrics}




