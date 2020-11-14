import sys
import model
import numpy as np

printStep = lambda text, step, freq : print('\t %s : %d'%(text, step)) if step % freq == 0 else None

class ModelParam:  # Simulation setup + model parameter
    def __init__(self):
        self.randGen = lambda dim, comp: np.array(np.random.uniform(low=0, high=1, size=(dim, dim)) < comp).astype(int)
        self.L = 75
        self.S = np.zeros(shape=(self.L, self.L)).astype(int)
        self.I = np.zeros_like(self.S)
        self.R = np.zeros_like(self.S)
        self.alpha = 5  # (m) lattice scale parameter
        self.infLT = 100 # (day) infectious life-time
        self.rho = 0.01  # tree density
        self.ell = 25
        self.epiC = int(self.L/2)
        self.S = self.randGen(dim=self.L, comp=self.rho)
        self.I[self.epiC, self.epiC] = 1
        self.S[np.where(self.I)] = 0
        self.model = ['exp', 'gaussian'][1]


class Settings:
    def __init__(self):
        self.plot = False
        self.save = False
        self.verbose = False
        self.pltFreq = 99
        self.ext='.png'

def R0_ensTseries(ensembleRepeats):
    import model
    import matplotlib.pyplot as plt
    """Run R0 finding simulation."""
    betas = [0.200]
    # ensembleAverage = np.zeros(ModelParam.infLT)
    R0_ensembe = np.zeros(shape=(ensembleRepeats, ModelParam().infLT))
    for beta in betas:
        print('Beta : {}'.format(beta))
        ModelParam.beta = beta
        for i in range(ensembleRepeats):
            printStep(text='Repeat ', step=i, freq=100)
            modelParam = ModelParam()
            settings = Settings()
            R0_tseries = np.zeros(modelParam.infLT)
            R0_tseries = model.findR0(pc=modelParam, numR=R0_tseries, settings=settings)
            R0_ensembe[i][:] = R0_tseries

    R0_average = R0_ensembe.sum(axis=0) / (i+1)
    plt.plot(R0_average, label='av')
    plt.legend()
    plt.show()
    return R0_average


if __name__ == '__main__':
    R0_average = R0_ensTseries(ensembleRepeats=100)
    import matplotlib.pyplot as plt
    # plt.plot(R0_average)
    # plt.show()
    sys.exit()