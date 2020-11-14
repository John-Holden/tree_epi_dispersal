from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import pickle

def ensAverage_data():
    print("Running: ensAverage_data...")
    SIR_data = pickle.load(open('SIR_ENSEMBLE.dat', 'rb'))
    SIR_ensembleAv = [np.zeros(10**4),np.zeros(10**4),np.zeros(10**4)]
    max_t = 0
    ens_size = 0
    for i in SIR_data:
        if 'S' in i:
            SIR_ensembleAv[0][0: len(SIR_data[i])] += SIR_data[i]
            ens_size += 1
            if len(SIR_data[i]) > max_t:
                max_t = len(SIR_data[i])
        if 'I' in i:
            SIR_ensembleAv[1][0: len(SIR_data[i])] += SIR_data[i]
            if len(SIR_data[i]) > max_t:
                max_t = len(SIR_data[i])
        if 'R' in i:
            SIR_ensembleAv[2][0: len(SIR_data[i])] += SIR_data[i]
            if len(SIR_data[i]) > max_t:
                max_t = len(SIR_data[i])

    SIR_ensembleAv[0] = SIR_ensembleAv[0][:max_t]
    SIR_ensembleAv[1] = SIR_ensembleAv[1][:max_t]
    SIR_ensembleAv[2] = SIR_ensembleAv[2][:max_t]
    print('...Done, ensemeble size {}'.format(ens_size))
    return np.array(SIR_ensembleAv) / ens_size

SIR_ensAv = ensAverage_data()
data = SIR_ensAv[1]
S0 = SIR_ensAv[0][0]
I0 = SIR_ensAv[1][0]
R0 = SIR_ensAv[2][0]
t_days = SIR_ensAv[0].shape[0]

def sumsq(p):
    beta, gamma = p
    def SIR(t,y):
        S = y[0]
        I = y[1]
        return([-beta*S*I, beta*S*I-gamma*I, gamma*I])
    sol = solve_ivp(SIR,[0, t_days], [S0, 1, 0],t_eval=np.arange(0,t_days, 0.5))
    return(sum( (sol.y[1][::2]-data)**2) )

msol = minimize(sumsq,[0.001,1],method='Nelder-Mead')
beta, gamma = msol.x
print('\t', msol, '\n............\n')
print('\t Found beta = {}, gamma = {}'.format(round(beta, 5), round(gamma, 5)))
def SIR(t,y):
    S = y[0]
    I = y[1]
    return([-beta*S*I, beta*S*I-gamma*I, gamma*I])

sol = solve_ivp(SIR,[0, t_days], [S0,I0,R0], t_eval=np.arange(0,t_days))
fig = plt.figure(figsize=(6,4))
plt.plot(sol.y[0],"b-")
plt.plot(sol.y[1],"r-")
plt.plot(sol.y[2],"g-")
plt.plot(np.arange(0, t_days, 10), SIR_ensAv[0][::10],"k*:", c='b', alpha=0.25)
plt.plot(np.arange(0, t_days, 10), SIR_ensAv[1][::10],"k*:", c='r', alpha=0.25)
plt.plot(np.arange(0, t_days, 10), SIR_ensAv[2][::10],"k*:", c='g', alpha=0.25)
plt.legend(["Susceptible","Infected","Reoved","Original Data"])
plt.savefig('EX_SIR_FITTING.pdf')
plt.show()
