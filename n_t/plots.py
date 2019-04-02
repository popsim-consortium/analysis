"""
Code for generating plots.
"""
<<<<<<< HEAD

import pandas
import seaborn as sns
=======
import pandas
>>>>>>> 692e85d433139f47d64784fc3736fa5f38f41ab3
import matplotlib
from matplotlib import pyplot as plt
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

sns.set_style("darkgrid")


def plot_stairway_Ne_estimate(infile, outfile):
    """
    figure of N(t) for single run of stairwayplot
    """
    nt = pandas.read_csv(infile, sep="\t", skiprows=5)
    nt = nt[nt['year'] > 10]
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    ax.plot(nt['year'], nt['Ne_median'], c="red")
    ax.plot(nt['year'], nt['Ne_2.5%'], c='grey')
    ax.plot(nt['year'], nt['Ne_97.5%'], c='grey')
    f.savefig(outfile, bbox_inches='tight')


def plot_compound_Ne_estimate(infiles, outfile):
    """
    figure of N(t) for multiple runs of stairwayplot
    """
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        nt = nt[nt['year'] > 10]
        ax.plot(nt['year'], nt['Ne_median'], c="red")
    f.savefig(outfile, bbox_inches='tight')


def plot_compound_smcpp(infiles, outfile):
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        ax.plot(nt['x'], nt['y'], c="red")
    f.savefig(outfile, bbox_inches='tight')


def plot_all_ne_estimates(sp_infiles, smcpp_infiles, outfile):
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    # plot smcpp estimates
    for infile in smcpp_infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line1, = ax.plot(nt['x'], nt['y'], c="red", alpha=0.8, label='smc++')
    # plot stairwayplot estimates
    for infile in sp_infiles:
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax.plot(nt['year'], nt['Ne_median'], c="blue", label='stairwayplot')
    # TODO add a plot  of true history

    ax.set_xlabel("time (years)")
    ax.set_ylabel("population size")
    ax.legend((line1, line2), ('smc++', 'stairwayplot'))
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)
