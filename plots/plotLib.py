import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
import numpy as np

pltParams = {'figure.figsize': (7.5, 5.5),
             'axes.labelsize': 15,
             'ytick.labelsize': 20,
             'xtick.labelsize': 20,
             'legend.fontsize': 'x-large'}
plt.rcParams.update(pltParams)

def frameLabel(step):
    if step < 10:
        return '000' + str(step)
    if step < 100:
        return '00' + str(step)
    if step < 1000:
        return '0' + str(step)
    if step < 10000:
        return str(step)

def pltSIR(S, I, R, dt):
    t = np.arange(0, len(S)) * dt
    plt.plot(t, S, c='green', label='S')
    plt.plot(t, I, c='red', label='I')
    plt.plot(t, R, c='black', label='R')
    plt.legend()
    plt.show()

def pltR0(R0_v_gen, save=False, print_out=True):
    if len(R0_v_gen) == 0:
        plt.plot([0, 1],  [0, 0], c='r', ls='--', label='Below threshold')
    else:
        t = np.arange(1, len(R0_v_gen)+1, 1)
        plt.plot(t, R0_v_gen)
        plt.scatter(t, R0_v_gen)
        plt.plot([1, t[-1]], [1, 1], c='r', ls='--')
    plt.ylabel(r'$R^i_0$')
    plt.xlabel('generation')
    if save:
        plt.savefig('r0_gen.pdf')
    plt.show()
    if print_out:
        p_out = {f'GEN {i}': R0_v_gen[i] for i in range(len(R0_v_gen))}
        print('\n RO(generation): \n', p_out)


def pltLG(S, I, R, dt):
    t = np.arange(0, len(S)) * dt

    plt.plot(t, S, c='green', label='S')
    plt.plot(t, I+R, c='red', label='I+R')
    plt.legend()
    plt.show()

def pltMaxD(maxD, dt):
    t = np.arange(0, len(maxD)) * dt
    plt.plot(t, maxD)
    plt.show()

def pltSim(S, I, R, t, anim, show):  # plot simulation time-steps
    pixSz = 15
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(np.where(I)[1], np.where(I)[0], s=pixSz, c='red')
    ax.scatter(np.where(S)[1], np.where(S)[0], s=pixSz, c='green')
    ax.scatter(np.where(R)[1], np.where(R)[0], s=pixSz, c='lightgray')
    ax.set_title( r'$t \approx ${} Days'.format(round(t, 3)), size=20)
    if anim:
        plt.savefig('./anim_dat/temp_frames/%s'%(frameLabel(step=t)))
    if show:
        plt.show()
    elif not show:
        plt.close()




