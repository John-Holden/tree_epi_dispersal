from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import sys

data = [1, 3, 6, 25, 73, 222, 294, 258, 237, 191, 125, 69, 27, 11, 4]

def SIR(t,y):
    S = y[0]
    I = y[1]
    return([-beta*S*I, beta*S*I-gamma*I, gamma*I])

def sumsq(p, SIR):
    beta, gamma = p
    def SIR(t,y):
        S = y[0]
        I = y[1]
        return([-beta*S*I, beta*S*I-gamma*I, gamma*I])
    sol = solve_ivp(SIR,[0,14],[762, 1, 0],t_eval=np.arange(0,14.2,0.2))
    print(sol.y[1][::5].shape)
    return(sum( (sol.y[1][::5]-data)**2) )

msol = minimize(sumsq,[0.001,1],method='Nelder-Mead')
beta,gamma = msol.x
sol = solve_ivp(SIR,[0, 14], [762,1,0],t_eval=np.arange(0,14.2,0.2))