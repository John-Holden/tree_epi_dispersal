import os
import math
import warnings
import numpy as np


from parameters_and_settings import ModelParamSet, Settings, Metrics


def time_print(time_seconds: int, msg: str = 'Simulation done in: ') -> str:
    """
    Pretty formatter to display time
    """
    seconds = math.floor(time_seconds)
    hrs = math.floor(seconds / 3600)
    mns = math.floor(seconds / 60)
    secs = seconds % 60

    if seconds < 60:
        print(f'{msg} {seconds} (s)')
    elif 60 <= seconds < 3600:
        print(f'{msg} {mns} (mins): {seconds} (s)')
    elif seconds >= 3600:
        print(f'{msg} {hrs} (Hrs): {mns%60} (mins): {secs} (s)')

    return f'{msg} {hrs} (Hrs): {mns%60} (mins): {secs} (s)'


def mk_new_dir(name: str, job_id:str):
    """
    Save new directory to file.
    """
    if job_id in ['1', 'local']:
        if os.path.exists(f'{os.getcwd()}/temp_dat_store/{name}') and 'HPC_MODE' in os.environ:
            raise FileExistsError(f'\n\t{os.getcwd()}/temp_dat_store/{name} \n'
                                  f'\tClear ensemble cache dumbass!!')

        elif os.path.exists(f'{os.getcwd()}/temp_dat_store/{name}') and not 'HPC_MODE' in os.environ:
            msg = f' \n\t Error, overwriting :{os.getcwd()}/temp_dat_store/{name} \n'
            warnings.warn(msg)
            return

        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}')
        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/info/')
        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/core_output')

    return


def save_meta_data(path_to_ensemble: str, job_id:str):
    """
    Save dispersal_model parameters, settings and metrics used in ensemble to file.
    """
    if job_id == '1' or job_id is None:  # write elapsed time to file
        np.save(f'{path_to_ensemble}/info/rhos', ModelParamSet.rhos)
        np.save(f'{path_to_ensemble}/info/betas', ModelParamSet.betas)
        with open(f'{path_to_ensemble}/info/ensemble_info.txt', 'w') as info_file:
            info_file.write("\n______Ensemble_______\n")
            info_file.write("\nNotes : '...' \n ")  # Note section to document results
            info_file.write("\n___Model parameters___\n")
            for param, value in vars(ModelParamSet).items():
                if param[0] == '_' or isinstance(value, staticmethod):
                    continue
                if param == 'alpha':
                    param, value = f'scale constant /{param}', f'{value} (m)'
                elif param == 'ell':
                    param, value = f'dispersal_type constant /{param}', f'{value} (m)'
                elif param == 'infLT':
                    param, value = f'infectious lifetime T', f'{value} (steps)'
                elif param == 'dispersal_model':
                    param = 'dispersal_type dispersal_model'
                elif param == 'r':
                    param = 'initial epicenter radius'
                elif param == 'init_n_infected':
                    param = 'initially infected'
                elif param == 'META_DATA' and value is None:
                    continue
                elif param == 'tend':
                    param, value = 'simulation run time', f'{value} (steps)'
                elif param == 'L':
                    try:
                        iter(value)
                        param = 'domain sizes L x L'
                        value = f'[{value} | len = {len(value)}]'
                    except TypeError:
                        param = 'domain size L x L'
                elif param == 'rho':
                    try:
                        iter(value)
                        param = 'tree densities /rho'
                        value = f'[{round(value[0], 3)}, {round(value[-1], 3)} | len = {len(value)}]'
                    except TypeError:
                        param = 'tree density /rho'
                        value = round(value, 3)
                elif param == 'beta':
                    try:
                        iter(value)
                        param = 'infectivity /beta'
                        value = f'[{round(value[0], 3)}, {round(value[-1], 3)} | len = {len(value)}]'
                    except TypeError:
                        param = 'infectivity /beta'
                        value = round(value, 3)

                info_file.write(f'\t - {param} : {value}\n')

            info_file.write("\n___Settings___\n")
            for param, value in vars(Settings).items():
                if param[0] == '_' or isinstance(value, staticmethod):
                    continue
                info_file.write(f'\t - {param} : {value}\n')

            info_file.write("\n___Metrics___\n")
            for param, value in vars(Metrics).items():
                if param[0] == '_' or isinstance(value, staticmethod):
                    continue
                info_file.write(f'\t - {param} : {value}\n')

            info_file.write("\n______Time elapsed_______\n")


def write_time(path_to_ensemble: str, time_elapsed: str, job_id:str):
    # write elapsed time
    with open(f'{path_to_ensemble}/info/ensemble_info.txt', 'a') as info_file:
        if job_id is None:
            info_file.write(f"\n {time_elapsed}\n")
            return

        info_file.write(f"core {job_id} | {time_elapsed}\n")

