"""
Code for generating plots.
"""
import pandas
import seaborn as sns
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


def plot_compound_msmc(infiles, outfile):
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        ax.plot(nt['x'], nt['y'], c="red")
    f.savefig(outfile, bbox_inches='tight')


def plot_all_ne_estimates(sp_infiles, smcpp_infiles, msmc_infiles, outfile):
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
    # plot msmc estimates
    for infile in msmc_infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line3, = ax.plot(nt['x'], nt['y'], c="purple", alpha=0.8, label='msmc')
    # TODO add a plot  of true history

    ax.set_xlabel("time (years)")
    ax.set_ylabel("population size")
    ax.legend((line1, line2, line3), ('smc++', 'stairwayplot', 'msmc'))
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)
