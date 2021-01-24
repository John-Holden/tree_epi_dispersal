import numpy as np
from timeit import default_timer as timer


class ModelParamSet:  # Set model parameter and simulation fields
    def __init__(self, rho, beta, L=1000, alpha=5):
        self.alpha = alpha  # (m) lattice scale parameter
        self.L = L  # L x L = domain dim : modelled size = alpha^2 * L^2 (m^2)
        self.infLT = 10 # (day) infectious life-time
        self.tend = 1000 # (day) final end time
        self.beta = beta # (day^{-1})) infectivity parameter
        self.rho = rho  # tree density
        self.ell = 1000  # (m) dispersal
        self.epiC = int(self.L/2)
        self.init_n_infected = 1 # control randomly distributed epicenters
        self.r = 0 # control centralised epicenter radius
        self.model = ['exp', 'gaussian'][1]


class Settings:  # Simulation setup
    def __init__(self):
        self.plot = True
        self.show = True
        self.anim = True
        self.pltFreq = 10
        self.verbose = 0
        self.ext='.png'
        self.boundary = False
        # if no `max_gen' infected trees left, end simulation
        self.gen_limit = 5 # R_0^i

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

if __name__ == '__main__':
    from runner_methods import run_lcl_R0_sensitivity, singleSim, R0_analysis
    # R0_analysis(singleSim(rho=0.01, beta=0.0001, L=1000)[1])
    run_lcl_R0_sensitivity(repeats=1, rho=0.01, beta=0.0001, alpha=7.5, box_sizes=[i for i in range(250, 2250, 250)])







