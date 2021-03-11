"""
Plot and analyse ensemble data.
"""
from parameters_and_settings import PATH_TO_DATA_STORE
from tree_epi_dispersal.plot_methods import plot_rho_beta_ensemble_1D
from tree_epi_dispersal.ensemble_analysis import process_R0_ensemble


if __name__ == '__main__':
    ens_name = f'{PATH_TO_DATA_STORE}2021-03-10-hpc-R0-vs-rho'
    ensemble, rhos, betas = process_R0_ensemble(ens_name, produce_landscape_control_package=True,
                                                process_nth_gen=3)

    plot_rho_beta_ensemble_1D(ensemble, rhos, betas)
