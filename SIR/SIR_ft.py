import numpy as np
import matplotlib.pyplot as plt
import sys


def field_upDate(S, I, beta, gamma, N):
    dS = - beta*S*I/N
    dI = beta*S*I/N - gamma*I
    dR = gamma*I
    return dS, dI, dR


def runSIR():
    # Return S,I,R simulations
    # Units hours
    N = 401
    I0, R0 = 2.0, 0.0
    S0 = N - I0 - R0
    # beta = 10 / (40 * 8 * 24)
    # gamma = 3 / (15 * 24)
    beta = 0.001
    gamma = 0.01
    tend = 25 # days
    dt = 0.1 # 6 mins
    steps  = int(tend*24/dt)
    S = np.zeros(steps)
    I = np.zeros_like(S)
    R = np.zeros_like(S)
    S[0], I[0], R[0] = S0, I0, R0
    for t in range(len(S)-1):
        S[t+1] = S[t] - dt*beta*S[t]*I[t]
        I[t+1] = I[t] + dt*beta*S[t]*I[t] - dt*gamma*I[t]
        R[t+1] = R[t] + dt*gamma*I[t]
        if t % 100 == 0:
            print(round(t * dt, 3), ' (hours) ', 'S =', round(S[t], 4))
    return S, I, R, np.linspace(0, tend, steps)


if __name__ == "__main__":
    S, I, R, t = runSIR()
    plt.plot(t, S, label='S')
    plt.plot(t, I, label='I')
    plt.plot(t, R, label='R')
    plt.ylabel('N')
    plt.xlabel('Days')
    plt.legend()
    plt.show()
    sys.exit()
