"""
Code for generating plots.
"""
import pandas
import seaborn as sns
import matplotlib
import msprime
import os
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
                            model, n_samp, generation_time, pop_id = 0, steps=None):

    #mass_migration_events = []
    #for demo in self.demographic_events:
    #if type(demo) == MassMigration:
    #    mass_migration_objects.append(demo)
    #    mass_migration_times.append(demo.time)
    
    ddb = msprime.DemographyDebugger(**model.asdict())
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1,end_time,end_time+1)
    pop_size = ddb.population_size_trajectory(steps=steps)
    # if mm abouve, switch pops at relavent time.
    pop_size = pop_size[:, pop_id]
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
        num_samples=num_samples, double_step_validation=False)
    steps = steps * generation_time

    f, ax = plt.subplots(1,3,sharex=True,sharey=True,figsize=(14, 7))
    # plot smcpp estimates
    for infile in smcpp_infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line1, = ax[0].plot(nt['x'], nt['y'], alpha=0.6)
    ax[0].plot(steps, 1/(2*coal_rate), c="black")
    ax[0].plot(steps, pop_size, c="#43464B")
    ax[0].set_title("smc++")

    # plot stairwayplot estimates
    for infile in sp_infiles:
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax[1].plot(nt['year'], nt['Ne_median'],alpha=0.6)
    ax[1].plot(steps, 1/(2*coal_rate), c="black")
    ax[1].plot(steps, pop_size, c="#43464B")
    ax[1].set_title("stairwayplot")

    # plot msmc estimates
    #msmc_sample_sizes = []
    for infile in msmc_infiles:
        fn = os.path.basename(infile)
        samp = fn.split(".")[0]
        #if samp not in msmc_sample_sizes:
        #    msmc_sample_sizes.append(samp)
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line3, = ax[2].plot(nt['x'], nt['y'], alpha=0.6)
    ax[2].plot(steps, 1/(2*coal_rate), c="black")
    ax[2].plot(steps, pop_size, c="#43464B")
    ax[2].set_title("msmc")

    for i in range(3):
        ax[i].set(xscale="log", yscale="log")
    ax[1].set_xlabel("time (years ago)")
    ax[0].set_ylabel("population size")
        #ax[i].legend((line1, line2, line3), ('smc++', 'stairwayplot', 'msmc'))
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)
