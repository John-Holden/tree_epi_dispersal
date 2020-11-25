import math


def timerPrint(seconds):
    dc = round(seconds - math.floor(seconds), 2)
    seconds = math.floor(seconds)
    hrs = math.floor(seconds / 3600)
    mns = math.floor(seconds / 60)
    secs = seconds % 60
    if seconds < 60:
        return '{} (s)'.format(round(secs + dc, 2))
    elif 60 <= seconds < 3600:
        return '{} (mins): {} (s)'.format(mns, secs)
    elif seconds >= 3600:
        return '{} (Hrs): {} (mins): {} (s)'.format(hrs, mns%60, secs)


def R0_history_sort(R0_trace):
    import numpy as np
    R0_count = np.zeros(10**4)
    num_trees_in_gen = np.zeros_like(R0_count)
    max_gen_in_sim = 0
    for site in R0_trace:
        inf_hist = R0_trace[site]
        R0_count[inf_hist[1]] += inf_hist[0]
        num_trees_in_gen[inf_hist[1]] += 1
        if inf_hist[1] > max_gen_in_sim:
            max_gen_in_sim = inf_hist[1]
    plt_maxGenn = max_gen_in_sim
    print('\n Max gen in sim {}'.format(max_gen_in_sim))
    R0_count = R0_count[:plt_maxGenn]
    num_trees_in_gen = num_trees_in_gen[:plt_maxGenn]
    mean_R0_for_gen = R0_count / num_trees_in_gen
    return mean_R0_for_gen