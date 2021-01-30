"""
Define classes for model parameters, settings and metrics.
"""
import numpy as np


class ModelParamSet:  # Set model parameter and simulation fields
    def __init__(self, rho, beta, L=1000, alpha=5, META_DATA=None):
        self.META_DATA = META_DATA
        self.alpha = alpha  # (m) lattice scale parameter
        self.L = L  # L x L = domain dim : modelled size = alpha^2 * L^2 (m^2)
        self.infLT = 50 # (day) infectious life-time
        self.tend = 100 # (day) final end time
        self.beta = beta # (day^{-1})) infectivity parameter
        self.rho = rho  # tree density
        self.ell = 1000  # (m) dispersal
        self.epiC = int(self.L/2)
        self.init_n_infected = 1 # control randomly distributed epicenters
        self.r = 0 # control centralised epicenter radius
        self.model = ['exp', 'gaussian'][1]


class Metrics:
    def __init__(self):
        self.endT = 0
        self.percT = None
        self.max_gen = None
        self.extinction = False
        self.extinctionT = None
        self.percolation = False
        self.mortality_ratio = 0
        self.R0 = np.zeros(3000)
        self.R0_histories = {}
        self.maxD = np.zeros(3000)
        self.numI = np.zeros(3000)
        self.numS = np.zeros(3000)
        self.numR = np.zeros(3000)


class Settings:  # Simulation setup
    plot = False
    show = False
    anim = False
    pltFreq = None
    verbose = 0  # verbosity
    ext = '.png'
    boundary = False  # Terminate upon infection reaching the boundary
    gen_limit = 1  # if no infected trees of order `gen-limit' are  left, end the simulation.





