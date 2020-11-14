import model
import numpy as np

class ModelParam:  # Simulation setup + model parameters
    # randGen: random numbers, compare to a value and return array type int
    def __init__(self, rho, beta):
        self.randGen = lambda L, comp: np.array(np.random.uniform(low=0, high=1, size=(L, L)) < comp).astype(int)
        self.L = 400
        self.S = np.zeros(shape=(self.L, self.L)).astype(int)
        self.I = np.zeros_like(self.S)
        self.R = np.zeros_like(self.S)
        self.alpha = 5  # (m) lattice scale parameter
        self.infLT = 100 # (day) infectious life-time
        self.tend = 31 * 4 # (day) final end time
        self.beta = beta # (day^{-1})) infectivity parameter
        self.rho = rho  # tree density
        self.ell = 25  # (m) dispersal co``````````nstant
        self.epiC = int(self.L/2)
        self.S = self.randGen(self.L, self.rho)
        self.r = 0 # epicenter radius
        if self.r == 0:
            self.I[self.epiC,self.epiC] = 1
        else:
            self.I[self.epiC-self.r:self.epiC+self.r,self.epiC-self.r:self.epiC+self.r] = 2
        self.S[np.where(self.I)] = 0
        self.model = ['exp', 'gaussian'][1]

class Settings:  # Simulation setup
    def __init__(self):
        self.plot = False
        self.show = True
        self.anim = True
        self.pltFreq = 31
        self.verbose1 = False
        self.verbose2 = False
        self.ext='.png'
        self.boundary = True

class Metrics:
    def __init__(self):
        self.endT = 0
        self.percolation = 0
        self.maxD = np.zeros(3000)
        self.numI = np.zeros(3000)
        self.numS = np.zeros(3000)
        self.numR = np.zeros(3000)

def singleSim(rho, beta):
    """
    Run a single instance of the model.
    """
    from plots.plotLib import pltSIR
    print('Running @ singleSim...')
    [modelParam, metrics] = model.runSim(pc=ModelParam(rho, beta), metrics=Metrics(), settings=Settings())
    pltSIR(S=metrics.numS, I=metrics.numI, R=metrics.numR, dt=1)
    print('....@ singleSim: DONE')
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
            vel_rpts, perc_rpts = run_ensemble(rho, beta, runs=N, MPrm=ModelParam, Mts=Metrics, Sts=Settings)
            vel_results[i, j] = vel_rpts
            perc_results[i, j] = perc_rpts

    saveFunc(vel_results, perc_results, rhos, betas, ModelParam, Settings(), ensName, jobId)
    print('....@GetPspace: DONE')
    return "Success"

def run_lcl_ens():
    import datetime
    from ensemble_methods import runVel_ensemble, mk_new_dir
    rhos = np.linspace(0.0, 0.10, 10)
    betas = [0.020, 0.03]
    repeats = 2
    date = datetime.datetime.today().strftime('%Y-%m-%d')
    ens_name = date + '-lcl-vel-ens'
    ens_name = mk_new_dir(ens_name)
    getPspace(runVel_ensemble, repeats, rhos, betas, ens_name, jobId='')


if __name__ == '__main__':
    import sys
    # run_lcl_ens()
    singleSim(rho=0.1, beta=0.1)
    sys.exit()


