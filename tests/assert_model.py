import unittest


from tree_epi_dispersal.model_dynamics_helpers import assert_correct_dispersal
from tree_epi_dispersal.model_dynamics import run_simulation
from tree_epi_dispersal.plot_methods import plot_test_dispersal

from parameters_and_settings import ParamsAndSetup


# class ExpectedGaussianl(unittest.TestCase):
    # def test_gaussian(self):
    #     print('TESTING GAUSSIAN')
    #
    #     test_scenario = {'rho':0.10, 'beta':0.01, 'ell':195, 'L':250, 'model':'gaussian'}
    #
    #     test_scenario['epi_c'] = int(test_scenario['L']/2)
    #     test_key = f"{test_scenario['epi_c']}{test_scenario['epi_c']}"
    #
    #     ParamsAndSetup['params'].tend = 1
    #     ParamsAndSetup['params'].L = test_scenario['L']
    #     ParamsAndSetup['params'].epi_center = test_scenario['epi_c']
    #
    #     ParamsAndSetup['params'].model = test_scenario['model']
    #     ParamsAndSetup['params'].ell = test_scenario['ell']
    #     assert_correct_dispersal()
    #
    #     ParamsAndSetup['setup'].plot = False
    #     ParamsAndSetup['setup'].verb = 0
    #
    #     ParamsAndSetup['metrics'].save_end_time = False
    #     ParamsAndSetup['metrics'].save_R0_history = True
    #     ParamsAndSetup['metrics'].save_time_series = False
    #     ParamsAndSetup['metrics'].save_mortality_ratio = False
    #
    #     observed_dispersal = []
    #     for i in range(1):
    #         sim_result = run_simulation(rho=test_scenario['rho'], beta=test_scenario['beta'], ell=test_scenario['ell'],
    #                                     test_mode=True)
    #
    #         R0_hist = sim_result['R0_hist']
    #         if i %10 == 0:
    #             print('Number infected : ', R0_hist[test_key][0])
    #         observed_dispersal.extend(R0_hist[test_key][2])
    #     plot_test_dispersal(test_scenario, observed_dispersal)
    #
    #     self.assertEqual(True, False)

class ExpectedGaussianl(unittest.TestCase):
    def test_power_law(self):
        print('TESTING POWER-LAW')
        test_scenario = {'rho':0.10, 'beta':1, 'ell':(200, 3.3), 'L':250, 'model':'power_law'}

        test_scenario['epi_c'] = int(test_scenario['L']/2)
        test_key = f"{test_scenario['epi_c']}{test_scenario['epi_c']}"

        ParamsAndSetup['params'].tend = 10
        ParamsAndSetup['params'].L = test_scenario['L']
        ParamsAndSetup['params'].epi_center = test_scenario['epi_c']

        ParamsAndSetup['params'].model = test_scenario['model']
        ParamsAndSetup['params'].ell = test_scenario['ell']
        assert_correct_dispersal()

        ParamsAndSetup['setup'].plot = False
        ParamsAndSetup['setup'].verb = 0

        ParamsAndSetup['metrics'].save_end_time = False
        ParamsAndSetup['metrics'].save_R0_history = True
        ParamsAndSetup['metrics'].save_time_series = False
        ParamsAndSetup['metrics'].save_mortality_ratio = False

        observed_dispersal = []
        for i in range(10):
            sim_result = run_simulation(rho=test_scenario['rho'], beta=test_scenario['beta'], ell=test_scenario['ell'],
                                        test_mode=True)

            R0_hist = sim_result['R0_hist']
            if i %10 == 0:
                print('Number infected : ', R0_hist[test_key][0])
            observed_dispersal.extend(R0_hist[test_key][2])
        plot_test_dispersal(test_scenario, observed_dispersal)

        self.assertEqual(True, False)

if __name__ == '__main__':
    unittest.main()
