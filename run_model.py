import model
import numpy as np
import datetime
from timeit import default_timer as timer
import sys


class ModelParamSet:  # Set model parameter and simulation fields
    def __init__(self, rho, beta):
        self.alpha = 5  # (m) lattice scale parameter
        self.L = 1000  # L = domain dim :: modelled size = (alpha * L)^2 (m^2)
        self.S = []
        self.I = []
        self.R = []
        self.infLT = 10 # (day) infectious life-time
        self.tend = 1000 # (day) final end time
        self.beta = beta # (day^{-1})) infectivity parameter
        self.rho = rho  # tree density
        self.ell = 1000  # (m) dispersal
        self.epiC = int(self.L/2)
        self.r = 0 # epicenter radius
        self.model = ['exp', 'gaussian'][1]
        self.randGen = lambda L, thresh: np.array(np.random.uniform(low=0, high=1, size=(L, L)) < thresh).astype(int)
        self.setField(epiLoc=['centralised', 'distributed'][1])  # Init SIR fields & epicenter

    def setField(self, epiLoc):
        self.S = self.randGen(self.L, thresh=self.rho)
        self.I = np.zeros_like(self.S)
        self.R = np.zeros_like(self.S)
        if epiLoc == 'centralised':
            self.setEpi()
        elif epiLoc == 'distributed':
            self.setRandEpi()

    def setEpi(self):
        if self.r == 0:
            self.I[self.epiC, self.epiC] = 1
        else:
            self.I[self.epiC-self.r:self.epiC+self.r, self.epiC-self.r:self.epiC+self.r] = 2
        self.S[np.where(self.I)] = 0

    def setRandEpi(self, n=10):
        randRow = np.array(np.random.randint(1, self.I.shape[0], size=n))
        randCol = np.array(np.random.randint(1, self.I.shape[0], size=n))
        randInf = tuple([randRow, randCol])
        self.I[randInf] = 2
        self.S[randInf] = 0

class Settings:  # Simulation setup
    def __init__(self):
        self.plot = False
        self.show = True
        self.anim = True
        self.pltFreq = 50
        self.verbose = False
        self.ext='.png'
        self.boundary = False

class Metrics:
    def __init__(self):
        self.endT = 0
        self.percT = None
        self.extinction = False
        self.extinctionT = None
        self.percolation = False
        self.mortality_ratio = 0
        self.R0 = np.zeros(3000)
        self.R0_trace = {}
        self.maxD = np.zeros(3000)
        self.numI = np.zeros(3000)
        self.numS = np.zeros(3000)
        self.numR = np.zeros(3000)

def singleSim(rho, beta):
    """
    Run a single instance of the model.
    """
    from helper_functions import timerPrint, R0_generation_mean
    from plots.plotLib import pltR0
    start = timer()
    print('Running @ singleSim...')
    [modelParam, metrics] = model.runSim(pc=ModelParamSet(rho, beta), metrics=Metrics(), settings=Settings())
    meanR0_vs_gen = R0_generation_mean(metrics.R0_trace)
    print('\n...Time elapsed = {}'.format(metrics.endT))
    print('...Percolation = {} @ time = {}'.format(metrics.percolation, metrics.percT))
    print('...Extinction = {} @ time = {}'.format(metrics.extinction, metrics.extinctionT))
    print('...Mortality Ratio = {} '.format(metrics.mortality_ratio))
    print('...Tot Number Removed = {} '.format(metrics.numR[-1]))
    elapsed = round(timer() - start, 2)
    print('\n@ singleSim DONE | {}'.format(timerPrint(elapsed)))
    pltR0(meanR0_vs_gen)
    return "Success"

def Pspace_iterator(run_ensemble, N, rhos, betas, ensName=None, jobId=None) -> "Success":
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations
    :param run_ensemble: method, type of ensemble running e.g. velocity/percolation/time-series
    :param rhos: float, av tree densities
    :param betas: float, infectivity constant
    """
    from helper_functions import timerPrint
    from ensemble_methods import save_ensemble, save_info

    start = timer()  # Time ensemble averaging process
    if jobId == '1' or jobId == 'local_test':
        save_info(ens_field_names=['R0_trace', 'extinction_time', 'mortality_ratio'],
                  rhos=rhos, betas=betas, param_set=ModelParamSet(0, 0), settings=Settings(),
                  ensemble_name=ensName, per_core_repeats=N)

    R0_results = np.zeros(shape=(len(betas), len(rhos), N))
    extinctionT_results = np.zeros(shape=(len(betas), len(rhos), N))
    mortality_ratio_results = np.zeros(shape=(len(betas), len(rhos), N))
    print('Running @Get Parameter Space...')
    for i, beta in enumerate(betas):
        print('i, ', i)
        for j, rho in enumerate(rhos):
            all_ens_fields = run_ensemble(rho, beta, runs=N, MPrm=ModelParamSet, Mts=Metrics, Sts=Settings)
            R0_results[i, j] = all_ens_fields[0]
            extinctionT_results[i, j] = all_ens_fields[1]
            mortality_ratio_results[i, j] = all_ens_fields[2]

    save_ensemble(ens_field_names=['R0_trace', 'extinction_time', 'mortality_ratio'],
                 ens_results=[R0_results, extinctionT_results, mortality_ratio_results],
                 ensemble_name=ensName, job_id=jobId)

    elapsed = round(timer() - start, 2)
    print('\n@ Ensemble Run DONE | {}'.format(timerPrint(elapsed)))
    return "Success"

def run_lcl_ens(repeats, rhos, betas):
    from ensemble_methods import runR0_ensemble, mk_new_dir
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-lcl-vel-ens'
    ens_name = mk_new_dir(ens_name)
    Pspace_iterator(runR0_ensemble, repeats, rhos, betas, ens_name, jobId='local_test')
    return 'Success'


if __name__ == '__main__':
    import sys
    run_lcl_ens(repeats=1, rhos=np.linspace(0.00, 0.01, 6), betas = np.linspace(0.001, 0.003, 2))
    # singleSim(rho=.05, beta=.003)
    sys.exit()


