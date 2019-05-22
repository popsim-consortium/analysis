"""
Code for generating plots.
"""
import pandas
import seaborn as sns
import matplotlib
import msprime
import os
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import numpy as np
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


def plot_all_ne_estimates(sp_infiles, smcpp_infiles, msmc_infiles, outfile,
                            model, n_samp, generation_time, species,
                            pop_id = 0, steps=None):

    ddb = msprime.DemographyDebugger(**model.asdict())
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1,end_time,end_time+1)
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
        num_samples=num_samples, double_step_validation=False)
    steps = steps * generation_time
    num_msmc = set([os.path.basename(infile).split(".")[0] for infile in msmc_infiles])
    num_msmc = sorted([int(x) for x in num_msmc])
    f, ax = plt.subplots(1,2+len(num_msmc),sharex=True,sharey=True,figsize=(14, 7))
    for infile in smcpp_infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line1, = ax[0].plot(nt['x'], nt['y'], alpha=0.8)
    ax[0].plot(steps, 1/(2*coal_rate), c="black")
    ax[0].set_title("smc++")
    for infile in sp_infiles:
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax[1].plot(nt['year'], nt['Ne_median'],alpha=0.8)
    ax[1].plot(steps, 1/(2*coal_rate), c="black")
    ax[1].set_title("stairwayplot")
    for i,sample_size in enumerate(num_msmc):
        for infile in msmc_infiles:
            fn = os.path.basename(infile)
            samp = fn.split(".")[0]
            if(int(samp) == sample_size):
                nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
                line3, = ax[2+i].plot(nt['x'], nt['y'],alpha=0.8)
        ax[2+i].plot(steps, 1/(2*coal_rate), c="black")
        ax[2+i].set_title(f"msmc, ({sample_size} samples)")
    plt.suptitle(f"{species}, population id {pop_id}", fontsize = 16)
    for i in range(2+len(num_msmc)):
        ax[i].set(xscale="log", yscale="log")
        ax[i].set_xlabel("time (years ago)")
    red_patch = mpatches.Patch(color='black', label='Coalescence rate derived Ne')
    ax[0].legend(frameon=False, fontsize=10, handles=[red_patch])
    ax[0].set_ylabel("population size")
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)
