"""
Define classes for dispersal_model parameters, settings and metrics.
"""
import os
import numpy as np
from tree_epi_dispersal.exceptions import InvalidDispersalSetup

PATH_TO_TEMP_STORE = f'{os.getcwd()}/temp_dat_store/'
PATH_TO_DATA_STORE = f'{os.getcwd()}/data_store/'


class ModelParamSet:  # Set dispersal_model parameter and simulation fields
    META_DATA = 'Exponentially distributed life-time dynamics : True'

    """ ________Default parameters________"""

    alpha = 5  # (m) lattice scale parameter
    L = 100  # L x L = domain dim : modelled size = alpha^2 * L^2 (m^2)
    infected_lt = 500  # (steps) infectious life-time -- exponentially-distributed mean = 1/T
    tend = 500  # (steps) final end time
    betas = [0.0005]  # (step^-1) infectivity-parameter
    rhos = [0.01]   # tree density
    # dispersal_type values from https://doi.org/10.1093/femsec/fiy049 are:
    #   - gaussian = 195m
    #   - power-law : a = 206, b=3.3
    # ell = (205, 3.3)  # (m) dispersal_type
    ell = 195
    epi_center = int(L/2)
    init_n_infected = 1  # control randomly distributed epicenters
    r = 0  # control centralised epicenter radius
    dispersal_model = ['exponential', 'gaussian', 'power_law'][1]
    model = ['SIR', 'ADB'][0]

    @staticmethod
    def update_epi_c():
        ModelParamSet.epi_center = int(ModelParamSet.L/2)

    @staticmethod
    def assert_correct_dispersal() -> bool:
        """ check dispersal_type configurations are correct """
        is_int_or_float = isinstance(ModelParamSet.ell, int) or isinstance(ModelParamSet.ell, float)
        is_gauss = ModelParamSet.dispersal_model == 'gaussian' and is_int_or_float
        is_exp = ModelParamSet == 'exponential' and is_int_or_float
        is_power_law = ModelParamSet.dispersal_model == 'power_law' and isinstance(ModelParamSet.ell, tuple) and len(
            ModelParamSet.ell) == 2

        valid = is_exp or is_gauss or is_power_law
        if not valid:
            raise InvalidDispersalSetup(ModelParamSet.dispersal_model, ModelParamSet.ell)

        return True

    @staticmethod
    def assert_config():
        assert int(2 * ModelParamSet.epi_center) == ModelParamSet.L, f'Incorrect epicenter {ModelParamSet.epi_center},' \
                                                                     f'L = {ModelParamSet.L}'
        assert ModelParamSet.assert_correct_dispersal()

    @staticmethod
    def adb_config():
        # set default ash dieback configuration
        ModelParamSet.fb_lt = [60, 21]  # mean and standard deviations for fruiting body life-time,
        ModelParamSet.I_tr_lt = 365 * 7
        ModelParamSet.model = 'ADB'


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
    plot_freq = 10
    verb = 2  # verbosity
    ext = '.png'
    percolation_bcd = False  # Terminate upon infection reaching the boundary
    max_generation_bcd = False  # if no infected trees of order @ max_gen are  left, end the simulation.

    @staticmethod
    def ensemble_config():
        # Initialise to set default attributes when running hpc simulation
        Settings.plot = False
        Settings.verb = 0


class EnsembleConfig:
    def __init__(self, config:str):
        if config == 'test':
            self.betas = [0.0001, 0.0005]
            self.rhos = [0.01, 0.02]

        elif config == 'full':
            rho_small = np.arange(0.0, 0.031, 0.0025)  # low-density trees
            rho_large = np.arange(0.04, 0.101, 0.01)
            self.rhos = np.hstack([rho_small, rho_large])  # len 20
            self.betas = np.arange(0.00000, 0.00032, 0.00002)  # len 16

        elif config == 'part':
            self.betas = [0.0001, 0.001]
            rho_small = np.arange(0.0, 0.031, 0.0025)  # low-density trees
            rho_large = np.arange(0.04, 0.101, 0.01)
            self.rhos = np.hstack([rho_small, rho_large])  # len 20

        else:
            raise TypeError


ParamsAndSetup = {'params': ModelParamSet, 'setup': Settings, 'metrics': Metrics, 'ensemble': EnsembleConfig}




