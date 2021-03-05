"""
Plot and analyse ensemble data.
"""
from parameters_and_settings import PATH_TO_DATA_STORE
from tree_epi_dispersal.plot_methods import plot_rho_beta_ensemble_1D
from tree_epi_dispersal.ensemble_analysis import process_R0_ensemble


if __name__ == '__main__':
    ens_name = f'{PATH_TO_DATA_STORE}/2021-01-30-hpc-R0-vs-rho'
    ensemble, rhos, betas = process_R0_ensemble(ens_name, field_of_interest='mean_R0_vs_gen_core_ensemble',
                                                produce_landscape_control_package=True)

    plot_rho_beta_ensemble_1D(ensemble, rhos, betas)