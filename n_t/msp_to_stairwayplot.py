import msprime as msp
import numpy as np
import allel
import argparse
import subprocess
from os.path import basename, splitext, join
# import sys

parser = argparse.ArgumentParser(
  description='Run stairwayplot on a tree sequence output from msprime.\
               Requires gnu-parallel and Stairwayplot including \
               the command-line version of Stairway_plot_output_summary,\
               which can be downloaded from: \
               https://sites.google.com/site/jpopgen/stairway-plot.')
parser.add_argument('infile',
                    help='path to tree sequence')
parser.add_argument('outdir',
                    help='directory for writing Stairwayplot input files')
parser.add_argument('nboots', type=int,
                    help='number of bootstrap replicates to run.')
parser.add_argument('stairwayplot_dir',
                    help="path to stairwayplot directory")
parser.add_argument('cores',
                    help="number of cores to use when fitting demographic \
                    models.")
parser.add_argument('seed', type=int,
                    help="random number seed for bootstrap resampling.")
parser.add_argument('mutation_rate', type=float,
                    help="mutation rate (per base per generation)")
parser.add_argument('generation_time', type=float,
                    help="generation time")
parser.add_argument('plot', type=bool,
                    help='plot Ne ~ t ? T/F')
args = parser.parse_args()

np.random.seed(args.seed)

# read in tree sequence
ts = msp.load(args.infile)

# count alleles, bootstrap over sites, return the SFS minus the 0% bin
haps = np.array(ts.genotype_matrix())
genotypes = allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
allele_counts = genotypes.count_alleles()
sfs = allel.sfs(allele_counts[:, 1])
sfs = sfs[1:len(sfs)]

# tmp directory for input files
command = ("cd " + args.outdir + ";" + "mkdir infiles")
subprocess.run(command, shell=True)

# write stairwayplot input
out = open((join(args.outdir, "infiles", splitext(basename(args.infile))[0]) +
            "_strwyplt.txt"), "w")
out.write(("msp" + "\t" + str(ts.num_samples) + "\t" +
           str(int(ts.sequence_length))+"\t" + str(1) + "\t" +
           str(ts.num_samples - 1) + "\n"))
# order is name,n_samples,sequence_length,lowest_sfs_bin,highest_sfs_bin
for x in sfs:
    out.write(str(int(x)) + "\t")
out.close()

# write bootstrapped inputs
for i in range(args.nboots):
    nsites = np.shape(allele_counts)[0]
    bootset = np.random.choice(np.arange(0, nsites, 1), nsites, replace=True)
    bootac = allele_counts[bootset, :]
    bootsfs = allel.sfs(bootac[:, 1])
    bootsfs = bootsfs[1:len(bootsfs)]
    out = open((join(args.outdir, "infiles",
                splitext(basename(args.infile))[0]) +
                "_strwyplt" + str(i) + ".txt"), "w")
    out.write(("msp" + "\t" + str(ts.num_samples) + "\t" +
               str(int(ts.sequence_length)) + "\t" +
               str(1) + "\t" + str(ts.num_samples-1) + "\n"))
    for x in bootsfs:
        out.write(str(int(x)) + "\t")
    out.close()

# fit models to bootstraps in parallel
command1 = ("cd " + args.stairwayplot_dir + ";" +
            "files=" + join(args.outdir, "infiles") + "/*;" +
            "parallel -j " + str(args.cores) + " " +
            "java -cp .:swarmops.jar " +
            "Stairway_plot_theta_estimation02 {1} 1 5000 " +
            "::: $files")
subprocess.run("echo $files", shell=True)
subprocess.run(command1, shell=True)

# moving output files...
command2 = ("cd " + args.outdir + ";" + "mkdir thetas;" + "mv " +
            join(args.outdir, "infiles/") +
            "*.addTheta thetas;" + "rm -r infiles")
subprocess.run(command2, shell=True)

# get median and 95% CI for Ne~t
command3 = ("cd " + args.stairwayplot_dir + ";" +
            "java Stairway_plot_output_summary_commandline "
            + join(args.outdir, "thetas") +
            " " + str(args.mutation_rate) + " " + str(args.generation_time) + " "
            + join(args.outdir, splitext(basename(args.infile))[1]
                   + "_estimated_Ne.txt") +
            ";" + "rm -r " + join(args.outdir, "thetas"))
subprocess.run(command3, shell=True)

if args.plot:
    from matplotlib import pyplot as plt
    import pandas
    nt = pandas.read_csv(join(args.outdir,
                              splitext(basename(args.infile))[0]+"_estimated_Ne.txt"),
                         sep="\t", skiprows=5)
    nt = nt[nt['year'] > 10]
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log",
           ylim=(np.min(nt['Ne_2.5%']-np.std(nt['Ne_2.5%'])),
                 np.max(nt['Ne_97.5%'])+np.std(nt['Ne_97.5%'])))
    ax.plot(nt['year'], nt['Ne_median'], c="red")
    ax.plot(nt['year'], nt['Ne_2.5%'], c='grey')
    ax.plot(nt['year'], nt['Ne_97.5%'], c='grey')
    f.savefig(join(args.outdir, splitext(basename(args.infile))[0]+"_Ne_est.png"),
              bbox_inches='tight')
