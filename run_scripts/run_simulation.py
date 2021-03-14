from tree_epi_dispersal.execute import single_sim
from parameters_and_settings import ParamsAndSetup

if __name__ == '__main__':
    ParamsAndSetup['params'].L = 100
    ParamsAndSetup['params'].update_epi_c()
    ParamsAndSetup['setup'].plot_freq = 1
    ParamsAndSetup['setup'].show = False

    single_sim(rho=0.02, beta=0.00050, ell=195, model='ADB')
