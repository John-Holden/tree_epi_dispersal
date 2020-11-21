import model
import numpy as np
import sys

class ModelParamSet:  # Set model parameter and simulation fields
    def __init__(self, rho, beta):
        self.alpha = 5  # (m) lattice scale parameter
        self.L = 1200  # L = domain dim :: modelled size = (alpha * L)^2 (m^2)
        self.S = []
        self.I = []
        self.R = []
        self.infLT = 10 # (day) infectious life-time
        self.tend = 200 # (day) final end time
        self.beta = beta # (day^{-1})) infectivity parameter
        self.rho = rho  # tree density
        self.ell = 1600  # (m) dispersal
        self.epiC = int(self.L/2)
        self.r = 1 # epicenter radius
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
            self.setRandEpi(n=10)

    def setEpi(self):
        if self.r == 0:
            self.I[self.epiC, self.epiC] = 1
        else:
            self.I[self.epiC-self.r:self.epiC+self.r, self.epiC-self.r:self.epiC+self.r] = 2
        self.S[np.where(self.I)] = 0

    def setRandEpi(self, n):
        randRow = np.array(np.random.randint(1, self.I.shape[0], size=n))
        randCol = np.array(np.random.randint(1, self.I.shape[0], size=n))
        randInf = tuple([randRow, randCol])
        self.I[randInf] = 2
        self.S[randInf] = 0

class Settings:  # Simulation setup
    def __init__(self):
        self.plot = True
        self.show = True
        self.anim = True
        self.pltFreq = 33
        self.verbose1 = False
        self.verbose2 = False
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
        self.maxD = np.zeros(3000)
        self.numI = np.zeros(3000)
        self.numS = np.zeros(3000)
        self.numR = np.zeros(3000)

def singleSim(rho, beta):
    """
    Run a single instance of the model.
    """
    from plots.plotLib import pltSIR, pltR0
    print('Running @ singleSim...')
    [modelParam, metrics] = model.runSim(pc=ModelParamSet(rho, beta), metrics=Metrics(), settings=Settings())
    print('...Mean R0 = {}'.format(round(metrics.R0.mean(), 3)))
    print('...Time elapsed = {}'.format(metrics.endT))
    print('...Percolation = {} @ time = {}'.format(metrics.percolation, metrics.percT))
    print('...Extinction = {} @ time = {}'.format(metrics.extinction, metrics.extinctionT))
    print('...Mortality Ratio = {} '.format(metrics.mortality_ratio))
    print('@ singleSim: DONE')
    # pltSIR(S=metrics.numS, I=metrics.numI, R=metrics.numR, dt=1)
    pltR0(R0=metrics.R0, perc=metrics.percolation, percT=metrics.percT, dt=1)
    return "Success"

def getPspace(run_ensemble, N, rhos, betas, ensName=None, jobId=None) -> "Success":
    """
    Get parameter space of model over rho/beta by ensemble-averaging simulations
    :param run_ensemble: method, type of ensemble running e.g. velocity/percolation/time-series
    :param rhos: float, av tree densities
    :param betas: float, infectivity constant
    """
    from ensemble_methods import saveFunc

    vel_results = np.zeros(shape=(len(betas), len(rhos), N))
    perc_results = np.zeros(shape=(len(betas), len(rhos), N))
    print('Running @GetPspace...')
    for i, beta in enumerate(betas):
        for j, rho in enumerate(rhos):
            vel_rpts, perc_rpts = run_ensemble(rho, beta, runs=N, MPrm=ModelParamSet, Mts=Metrics, Sts=Settings)
            vel_results[i, j] = vel_rpts
            perc_results[i, j] = perc_rpts

    saveFunc(vel_results, perc_results, rhos, betas, ModelParamSet, Settings(), ensName, jobId)
    print('....@GetPspace: DONE')
    return "Success"

def run_lcl_ens():
    import datetime
    from timeit import default_timer as timer
    from ensemble_methods import runVel_ensemble, mk_new_dir
    rhos = np.linspace(0.04, 0.04, 1)
    betas = np.linspace(0.3, 0.3, 1)
    repeats = 1
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-lcl-vel-ens'
    ens_name = mk_new_dir(ens_name)
    st = timer()
    getPspace(runVel_ensemble, repeats, rhos, betas, ens_name, jobId='')
    print('\t Took: {} (s)'.format(timer() - st))


if __name__ == '__main__':
    import sys
    # run_lcl_ens()
    for i in range(1):
        singleSim(rho=.01, beta=.00005)
    sys.exit()


