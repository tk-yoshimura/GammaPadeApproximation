import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dirpath_summary = '../results/'
suffix_summary = '.csv'

dirpath_figure = '../figures/'
suffix_figure = '.svg'

targets = [
    'gamma_approx',
    'loggamma_approx',
    'digamma_approx',
    'invgamma_approx',
    'polygamma1_approx',
    'polygamma2_approx',
    'polygamma3_approx',
    'polygamma4_approx',
    'polygamma5_approx',
    'polygamma6_approx',
    'polygamma7_approx',
    'polygamma8_approx',
    'polygamma9_approx',
    'polygamma10_approx',
    'polygamma11_approx',
    'polygamma12_approx',
    'polygamma13_approx',
    'polygamma14_approx',
    'polygamma15_approx',
    'polygamma16_approx',
]

options = [
    ['logy'],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
]

for target, option in zip(targets, options):
    data = pd.read_csv(dirpath_summary + target + suffix_summary) 

    x, y, err = data['x'], data['y_actual'], data['error(rate)']
    err = np.where(np.abs(y) < 1e-305, 0, err)

    plt.clf()
    fig = plt.figure(figsize=(12, 6))

    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

    ax1.plot(x, y)
    ax1.grid()
    ax1.set_ylabel('actual')

    ax2.plot(x, err)
    ax2.grid()
    ax2.set_ylabel('error(rate)')

    if 'logx' in option:
        ax1.set_xscale('log')
        ax2.set_xscale('log')

    if 'logy' in option:
        ax1.set_yscale('log')

    if target.startswith('polygamma'):
        ax1.set_ylim([max(np.min(y)-0.1, -5), min(np.max(y)+0.1, 5)])
        
    plt.savefig(dirpath_figure + target + suffix_figure, bbox_inches='tight')
