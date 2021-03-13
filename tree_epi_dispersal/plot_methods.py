import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Iterable, Any
from parameters_and_settings import PATH_TO_TEMP_STORE, ParamsAndSetup

pltParams = {'figure.figsize': (7.5, 5.5),
             'axes.labelsize': 15,
             'ytick.labelsize': 20,
             'xtick.labelsize': 20,
             'legend.fontsize': 'x-large'}
plt.rcParams.update(pltParams)


# -------------------Simulation plotting methods-------------------#
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
        plt.plot([0, 1], [0, 0], c='r', ls='--', label='Below threshold')
    else:
        t = np.arange(1, len(R0_v_gen) + 1, 1)
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
    plt.plot(t, I + R, c='red', label='I+R')
    plt.legend()
    plt.show()


def pltMaxD(maxD, dt):
    t = np.arange(0, len(maxD)) * dt
    plt.plot(t, maxD)
    plt.show()


def plt_sim_frame(S, I, R, t, save, show, msg=None):  # plot simulation time-steps
    pixSz = 8
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(np.where(I)[1], np.where(I)[0], s=pixSz, c='red')
    ax.scatter(np.where(S)[1], np.where(S)[0], s=pixSz, c='green')
    ax.scatter(np.where(R)[1], np.where(R)[0], s=pixSz, c='lightgray')

    title = r'$t \approx ${} Days '.format(round(t, 3))
    title = title + msg if msg else title
    ax.set_title(title, size=20)
    if save:  # save to animations folder
        plt.savefig('./anim_dat/temp_frames/%s' % (frameLabel(step=t)))
    if show:
        plt.show()
    elif not show:
        plt.close()


def plt_adb_frame(S: tuple, E: tuple, I_fb: list, R_fb: list, t: Any):

    pixSz = 12.5
    fig, ax = plt.subplots(figsize=(7, 7))
    plt.title(f'T : {t}')
    ax.scatter(I_fb[1], I_fb[0], s=pixSz, c='red')
    ax.scatter(R_fb[1], R_fb[0], s=pixSz, c='grey')
    ax.scatter(S[1], S[0], s=pixSz, c='green')
    ax.scatter(E[1], E[0], s=pixSz, c='orange')
    plt.show()
    plt.close()


# ------------------- Ensemble plotting -------------------#
def plot_rho_beta_ensemble_1D(ensemble: np.ndarray, rhos: np.ndarray, betas: np.ndarray):
    """
    For an ensemble, plot a series of 1D lines, rho axis, for multiple beta values.
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    for i in range(len(ensemble)):
        rho_line = ensemble[i]
        ax.plot(rhos, rho_line, label=f'beta = {round(betas[i], 5)}')
        ax.scatter(rhos, rho_line)
    print('BETAS: ', betas)
    ax.plot([0, rhos[-1]], [1, 1], color='r', ls='--')
    plt.legend()
    plt.show()


def plot1D_mean(rhos, betas, ens_mean, save=False):
    for i, beta in enumerate(betas):
        plt.plot(rhos, ens_mean[i], c='C' + str(i),
                 label=r'$\beta$ = {}'.format(round(beta, 8)))
        plt.scatter(rhos, ens_mean[i], c='C' + str(i), )

    plt.hlines(y=1, xmin=rhos[0], xmax=rhos[-1], ls='--')
    plt.legend()
    if save:
        plt.savefig('ens_dat.pdf')
    plt.show()


def plot_R0_ens_vs_gen(path_to_ensemble: str, ensemble_mean: dict, save=False):
    """
    Process ensemble_core_averaged R0 history and plot generational mean R0
    """
    from tree_epi_dispersal.model_dynamics_helpers import avg_multi_dim
    box_sizes = np.load(f'{path_to_ensemble}/info/box_sizes.npy')[::-1]
    c = 0
    for box_size in ensemble_mean:
        assert int(box_size) in box_sizes
        c += 1
    assert c == len(box_sizes)
    R0_v_L = np.zeros(len(box_sizes))

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for i, box_size in enumerate(box_sizes):
        R0_vs_gen = avg_multi_dim(ensemble_mean[str(box_size)])[:21]
        R0_v_L[i] = R0_vs_gen[0]
        gen = np.arange(1, len(R0_vs_gen) + 1)
        ax.plot(gen, R0_vs_gen, label=f'{box_size}  {box_size}')
        ax.scatter(gen, R0_vs_gen)

    ax.set_xlim(0.5, 10.5)
    plt.legend()
    if save:
        np.save('R0_vs_L_alpha_()m', R0_v_L)
        plt.tight_layout()
        plt.savefig('R0_vs_gen.pdf')
    plt.show()


def plot_R0_ens_vs_L(save=False):
    """
    For different alpha values, plot saved R0 values against domain size L
    """
    data_sets = {'2021-01-24-hpc-R0-generation-alpha-5': '5',
                 '2021-01-24-hpc-R0-generation-alpha-10': '10',
                 '2021-01-24-hpc-R0-generation-alpha-7_5': '7_5'}

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for data_set in data_sets:
        box_sizes = np.load(f'{PATH_TO_TEMP_STORE}/{data_set}/info/box_sizes.npy')[::-1]
        alpha_5 = np.load(f'{PATH_TO_TEMP_STORE}/{data_set}/R0_vs_L_alpha_{data_sets[data_set]}m.npy')
        ax.plot(box_sizes, alpha_5, label=rf'$\alpha = $ {data_sets[data_set].replace("_", ".")}')
        ax.scatter(box_sizes, alpha_5)

    plt.legend()
    plt.tight_layout()
    if save:
        plt.savefig('R0_vs_L_vs_alpha.pdf')
    plt.show()


def plot_test_dispersal(test_scenario: dict, actual_dispersal: Iterable):
    """
    Plot the expected kernel against the observed kernel.
    :param model:
    :param ell:
    :param actual_dispersal:
    :return:
    """
    import numpy as np
    import seaborn
    from tree_epi_dispersal.model_dynamics_helpers import model_selector

    # plt.hist(actual_dispersal)
    seaborn.displot(actual_dispersal)
    plt.ylabel('count')
    plt.xlabel('dist')
    plt.show()

    beta = test_scenario['beta']
    ell = test_scenario['ell']
    L = test_scenario['L']

    analytic_model = model_selector()
    dist = np.linspace(-L / 2, L / 2, 100) * ParamsAndSetup['params'].alpha
    center = test_scenario['epi_c'] * ParamsAndSetup['params'].alpha
    plt.plot(dist, analytic_model(dist, beta, ell))
    plt.ylabel('pr')
    plt.xlabel('dist')
    plt.show()
