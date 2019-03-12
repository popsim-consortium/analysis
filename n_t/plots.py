"""
Code for generating plots.
"""

import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas

def plot_stairway_Ne_estimate(infile, outfile):

    nt = pandas.read_csv(infile, sep="\t", skiprows=5)
    nt = nt[nt['year'] > 10]
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    # JK: turned off these limits as they weren't working for me.
           # ylim=(np.min(nt['Ne_2.5%']-np.std(nt['Ne_2.5%'])),
           #       np.max(nt['Ne_97.5%'])+np.std(nt['Ne_97.5%'])))
    ax.plot(nt['year'], nt['Ne_median'], c="red")
    ax.plot(nt['year'], nt['Ne_2.5%'], c='grey')
    ax.plot(nt['year'], nt['Ne_97.5%'], c='grey')
    f.savefig(outfile, bbox_inches='tight')
