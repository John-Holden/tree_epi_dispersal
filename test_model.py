import unittest
import numpy as np
from PARAMETERS_AND_SETUP import ModelParamSet, Metrics, Settings
from helper_methods import setFields

class TestStringMethods(unittest.TestCase):

    def test_expected_number_of_infecteds_upon_init(self):
        # test different ways initially infected trees can be setup
        number = None  # case 1, single centralised epicenter
        radius = 0
        S,I,R=setFields(L=1000, rho=0.01, epiC=500, r=radius, init_n_infected=number)
        assert len(np.where(I)[0]) == 1

        number = None
        radius = 2
        S, I, R = setFields(L=1000, rho=0.01, epiC=500, r=radius, init_n_infected=number)
        assert len(np.where(I)[0]) == 16

        number = 10
        radius = None
        S, I, R = setFields(L=1000, rho=0.01, epiC=500, r=radius, init_n_infected=number, epiType='distributed')
        assert len(np.where(I)[0]) == 10

    def expected_number_of_trees_over_domain(self):
        # does the number trees match the expected ?
        rho = 0.01
        L = 1000
        alpha = 5
        model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
        print(f'Modelled area = {L*alpha/1000}km x {L*alpha/1000}km')
        actual_tree_number = model_params.S.sum()
        expected_tree_number = rho*(L)**2
        print(f'av = {actual_tree_number}')
        print(f'expected_tree_number = {expected_tree_number}')
        assert actual_tree_number == expected_tree_number  # /Error expected same number of trees

    def expected_number_of_trees_for_different_alpha_values(self):
        from model_dynamics import ijDistance
        rho = 0.01
        L = 1000
        alpha = 5 # m
        model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
        print(f'alpha = {alpha}m')
        print(f'Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
        actual_tree_number_alpha_5m = model_params.S.sum()
        print(f'alpha=5m | actual tree number = {actual_tree_number_alpha_5m}')
        dist_arr = np.array([i for i in range(len(model_params.S[0]))]).astype(float)
        array_dist = ijDistance([0, 0], [0, dist_arr[-1]])
        expected_size_km = L * alpha / 1000
        actual_size_km = ijDistance([0, 0], [0, dist_arr[-1] + 1]) * alpha / 1000
        assert int(array_dist) == L - 1  # /Error expect correct domain size
        assert expected_size_km == actual_size_km  # /Error expect correct domain size
        dist_arr -= len(model_params.S[0]) / 2
        epi_c = dist_arr[model_params.epiC]
        assert epi_c == 0
        assert dist_arr[-1] == int(L / 2) - 1
        print(f'Domain size = {actual_size_km}km')

    def test_different_L_and_alpha_values_which_give_same_modelled_area(self):
        from helper_methods import setFields
        rho = 0.01
        L = 1000
        alpha = 5  # (m)
        model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
        print(f'alpha = {alpha}m, L = {L}, Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
        #-----------
        rho = 0.01
        L = 500
        alpha = 10  # (m)
        model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
        # setFields(L, rho, epiC, r, init_n_infected)
        actual_tree_number_alpha_5m = model_params.S.sum()
        expected_tree_number5m = rho * L ** 2
        print(f'alpha = {alpha}m, L = {L}, Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
        actual_tree_number_alpha_10m = model_params.S.sum()
        expected_tree_number10m = rho * L**2
        print(f'actual tree number alpha=5m: {actual_tree_number_alpha_5m}, '
              f'actual tree number alpha=10m: {actual_tree_number_alpha_10m}')
        print(f'expected tree number alpha=5m: {expected_tree_number5m}, '
              f'expected tree number alpha=10m: {expected_tree_number10m}')

if __name__ == '__main__':
    unittest.main()


