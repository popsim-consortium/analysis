"""
Code for generating plots.
"""
import pandas
import os
import matplotlib
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
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


def plot_compound_smcpp(infiles, outfile, model, n_samp, generation_time, pop_id=0):
    coal_rate, P, steps = popn_coal_rate(model, pop_id, n_samp, generation_time)
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        ax.plot(nt['x'], nt['y'], c="red")
    ax.plot(steps, 1/(2*coal_rate), c="black")
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
                          pop_id=0, steps=None):
    ddb = model.get_demography_debugger()
    if steps is None:
        if(len(ddb.epochs) == 1):
            end_time = 100000
        else:
            end_time = int(ddb.epochs[-2].end_time) + 10000
        steps = np.linspace(1, end_time, 1000)
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   num_samples=n_samp,
                                                   double_step_validation=False)
    mm = [x for x in ddb.demographic_events if x.type == "mass_migration"]
    census_size = ddb.population_size_trajectory(steps=steps)
    for m in reversed(mm):
        if m.proportion == 1.0:
            n = (steps > m.time)
            census_size[n, m.source] = census_size[n, m.dest]
        else:
            print("Error: census size estimate requires that MassMigration proportion == 1.0")
            sys.exit(1)
    steps = steps * generation_time
    num_msmc = set([os.path.basename(infile).split(".")[0] for infile in msmc_infiles])
    num_msmc = sorted([int(x) for x in num_msmc])
    f, ax = plt.subplots(1, 2+len(num_msmc), sharex=True, sharey=True, figsize=(14, 7))

    outLines = []
    for i, infile in enumerate(smcpp_infiles):
        nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
        line1, = ax[0].plot(nt['x'], nt['y'], alpha=0.8)
        for j in range(len(nt["x"])):
            outLines.append([nt["x"][j], nt["y"][j], "smcpp", "r" + str(i + 1)])
    ax[0].plot(steps, 1/(2*coal_rate), c="black")
    ax[0].set_title("smc++")

    for i, infile in enumerate(sp_infiles):
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax[1].plot(nt['year'], nt['Ne_median'], alpha=0.8)
        for j in range(0, len(nt["year"]), 2):
            outLines.append([nt["year"][j], nt["Ne_median"][j], "sp", "r" + str(i + 1)])

    ax[1].plot(steps, 1/(2*coal_rate), c="black")
    ax[1].set_title("stairwayplot")
    for i, sample_size in enumerate(num_msmc):
        ct = 0
        for infile in msmc_infiles:
            fn = os.path.basename(infile)
            samp = fn.split(".")[0]
            if(int(samp) == sample_size):
                nt = pandas.read_csv(infile, usecols=[1, 2], skiprows=0)
                line3, = ax[2+i].plot(nt['x'], nt['y'], alpha=0.8)
                for j in range(len(nt["x"])):
                    outLines.append([nt["x"][j], nt["y"][j], "msmc_" +
                                     str(sample_size), "r" + str(ct + 1)])
                ct += 1
        ax[2+i].plot(steps, 1/(2*coal_rate), c="black")
        ax[2+i].set_title(f"msmc, ({sample_size} samples)")
    plt.suptitle(f"{species}, population id {pop_id}", fontsize=16)
    for i in range(2+len(num_msmc)):
        ax[i].set(xscale="log")
        ax[i].set_xlabel("time (years ago)")
    red_patch = mpatches.Patch(color='black', label='Coalescence rate derived Ne')
    ax[0].legend(frameon=False, fontsize=10, handles=[red_patch])
    ax[0].set_ylabel("population size")
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)

    txtOUT = os.path.join(os.path.dirname(outfile),"_".join([species,model.id,"pop"+str(pop_id),"sizes"])+".txt")
    with open(txtOUT, "w") as fOUT:
        fOUT.write("\t".join([str(x) for x in ["x", "y", "method", "rep"]]) + "\n")
        for i in range(len(steps)):
            fOUT.write("\t".join([str(x) for x in [steps[i],
                                                   1/(2*coal_rate[i]),
                                                   "coal",
                                                   "r1"]]) + "\n")
            fOUT.write("\t".join([str(x) for x in [steps[i],
                                                   census_size[i][pop_id],
                                                   "census",
                                                   "r1"]]) + "\n")
        for i in range(len(outLines)):
            fOUT.write("\t".join([str(x) for x in outLines[i]])+"\n")


def plot_stairwayplot_coalrate(sp_infiles, outfile,
                               model, n_samp, generation_time, species,
                               pop_id=0, steps=None):  # JRA

    ddb = model.get_demography_debugger()
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1, end_time, end_time+1)
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    print(num_samples)

    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   num_samples=num_samples,
                                                   double_step_validation=False)
    steps = steps * generation_time
    f, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(7, 7))
    ax.plot(steps, 1/(2*coal_rate), c="black")
    for infile in sp_infiles:
        nt = pandas.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax.plot(nt['year'], nt['Ne_median'], alpha=0.8)
    ax.plot(steps, 1/(2*coal_rate), c="black")
    ax.set_title("stairwayplot")
    plt.suptitle(f"{species}, population id {pop_id}", fontsize=16)
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("time (years ago)")
    red_patch = mpatches.Patch(color='black', label='Coalescence rate derived Ne')
    ax.legend(frameon=False, fontsize=10, handles=[red_patch])
    ax.set_ylabel("population size")
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)


def popn_coal_rate(model, pop_id, n_samp, generation_time, steps=None):
    """
    returns tuple (coal_rate, P, steps) for pop_id
    conditional on the model and configuration
    """
    ddb = model.get_demography_debugger()
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1, end_time, end_time+1)
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   num_samples=num_samples,
                                                   double_step_validation=False)
    steps = steps * generation_time
    return coal_rate, P, steps


def plot_sfs(s, outfile):
    """
    Plot the SFS for this simulation
    """
    bins = [n + 1 for n in range(len(s[0]))]
    vals = []
    for i in range(len(s)):
        vals.append([int(x) for x in s[i]])
    if len(s) == 2:
        f, ax = plt.subplots(1, 2, sharey=True, tight_layout=True, figsize=(8, 3))
        ax[0].bar(bins, vals[0])
        ax[1].bar(bins, vals[1])
        ax[0].set_title("all sites")
        ax[1].set_title("masked sites")
        ax[0].set_ylabel("counts")
        ax[0].set_xlabel("derived allele frequency")
        ax[1].set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()
    else:
        f, ax = plt.subplots(figsize=(3, 3))
        ax.bar(bins, vals[0])
        ax.set_title("all sites")
        ax.set_ylabel("counts")
        ax.set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()
