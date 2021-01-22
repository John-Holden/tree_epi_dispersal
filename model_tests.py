import unittest
import numpy as np
from run_model import ModelParamSet, Metrics, Settings

def expected_number_of_trees_over_domain():
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

def expected_number_of_trees_for_different_alpha_values():
    from model import ijDistance
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
    #-------------
    rho = 0.01
    L = 1000
    alpha = 10 #m
    model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
    print(f'alpha = {alpha}m')
    print(f'Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
    actual_tree_number_alpha_10m = model_params.S.sum()
    assert actual_tree_number_alpha_5m == actual_tree_number_alpha_10m
    print(f'alpha=10m | actual tree number = {actual_tree_number_alpha_10m}')
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
    # -------------
    rho = 0.01
    L = 1000
    alpha = 20  # m
    model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
    print(f'alpha = {alpha}m')
    print(f'Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
    actual_tree_number_alpha_20m = model_params.S.sum()
    print(f'alpha=20m | actual tree number = {actual_tree_number_alpha_20m}')
     # /Error expected same number of trees
    assert actual_tree_number_alpha_20m == actual_tree_number_alpha_5m  # /Error expected same number of trees
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
    #-------------

def test_different_L_and_alpha_values_which_give_same_modelled_area():
    rho = 0.01
    L = 1000
    alpha = 5  # (m)
    model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
    print(f'alpha = {alpha}m, L = {L}, Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
    actual_tree_number_alpha_5m = model_params.S.sum()
    expected_tree_number5m = rho * L**2
    #-----------
    rho = 0.01
    L = 500
    alpha = 10  # (m)
    model_params = ModelParamSet(rho=rho, beta=0.001, L=L, alpha=alpha)
    print(f'alpha = {alpha}m, L = {L}, Modelled area = {L * alpha / 1000}km x {L * alpha / 1000}km')
    actual_tree_number_alpha_10m = model_params.S.sum()
    expected_tree_number10m = rho * L**2
    print(f'actual tree number alpha=5m: {actual_tree_number_alpha_5m}, '
          f'actual tree number alpha=10m: {actual_tree_number_alpha_10m}')
    print(f'expected tree number alpha=5m: {expected_tree_number5m}, '
          f'expected tree number alpha=10m: {expected_tree_number10m}')

if __name__ == '__main__':
    print('-----------------------------')
    print('Testing domain density is as expected')
    print('-----------------------------')
    expected_number_of_trees_over_domain()
    print('-----------------------------')
    print('Testing the domain over different alpha values, one L value')
    print('-----------------------------')
    expected_number_of_trees_for_different_alpha_values()
    print('-----------------------------')
    print('Testing the domain over different alpha and L value to give same domain')
    print('-----------------------------')
    test_different_L_and_alpha_values_which_give_same_modelled_area()
