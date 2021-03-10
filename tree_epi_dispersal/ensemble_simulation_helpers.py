import os
import math

from warnings import warn

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
    if job_id == '1' or job_id is None:
        dir1 = f'{os.getcwd()}/temp_dat_store/{name}'
        sub_dir1 = f'{os.getcwd()}/temp_dat_store/{name}/info/'
        sub_dir2 = f'{os.getcwd()}/temp_dat_store/{name}/core_output'

        if os.path.exists(dir1):
            msg = f'FileExistsWarn: {os.getcwd()}/temp_dat_store/{name}'
            # raise FileExistsError(F'{os.getcwd()}/temp_dat_store/{name}')
            warn(msg)
            assert os.path.exists(sub_dir1)
            assert os.path.exists(sub_dir2)
            return

        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}')
        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/info/')
        os.mkdir(f'{os.getcwd()}/temp_dat_store/{name}/core_output')

    return


def save_meta_data(path_to_ensemble: str, job_id:str):
    """
    Save model parameters, settings and metrics used in ensemble to file.
    """
    if job_id == '1' or job_id is None:  # write elapsed time to file
        with open(f'{path_to_ensemble}/info/ensemble_info.txt', 'w') as info_file:
            info_file.write("\n______Ensemble_______\n")
            info_file.write("\nNotes : '...' \n ")  # Note section to document results
            info_file.write("\n___Model parameters___\n")
            for param, value in vars(ModelParamSet).items():
                if param == 'alpha':
                    param, value = f'scale constant /{param}', f'{value} (m)'
                elif param == 'ell':
                    param, value = f'dispersal constant /{param}', f'{value} (m)'
                elif param == 'infLT':
                    param, value = f'infectious lifetime T', f'{value} (steps)'
                elif param == 'model':
                    param = 'dispersal model'
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
                if param[0] == '_':
                    continue
                info_file.write(f'\t - {param} : {value}\n')

            info_file.write("\n___Metrics___\n")
            for param, value in vars(Metrics).items():
                if param[0] == '_':
                    continue
                info_file.write(f'\t - {param} : {value}\n')

            info_file.write("\n______Time elapsed_______\n")


def write_time(path_to_ensemble: str, time_elapsed: str, job_id:str):
    # write elapsed time
    with open(f'{path_to_ensemble}/info/ensemble_info.txt', 'a') as info_file:
        if job_id is None:
            info_file.write(f"\n {time_elapsed}\n")
            return

        info_file.write(f"\n core {job_id} | {time_elapsed}\n")

