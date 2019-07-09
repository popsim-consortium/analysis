"""
Module to run the required simulations
and create unlinked-sites mask
"""

import numpy as np
import msprime as msp
import pickle
import threading
import math
import allel
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
sns.set_style("darkgrid")


def plot_sfs(s, outfile):
    """
    Plot the SFS for this simulation
    """
    bins = [n + 1 for n in range(len(s[0]))]
    vals = []
    for i in range(len(s)):
        vals.append([int(x) for x in s[i]])
    if len(s) == 2:
        f, ax = plt.subplots(1,2,sharey=True,tight_layout=True,figsize=(8,3))
        ax[0].bar(bins, vals[0])
        ax[1].bar(bins, vals[1])
        ax[0].set_title("all sites")
        ax[1].set_title("ld-filtered sites")
        ax[0].set_ylabel("counts")
        ax[0].set_xlabel("derived allele frequency")
        ax[1].set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
    else:
        f, ax = plt.subplots(figsize=(3,3))
        ax.bar(bins, vals[0])
        ax.set_title("all sites")
        ax.set_ylabel("counts")
        ax.set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')


def unlinked(ts,thresh,num_threads):
    """
    Returns a boolean mask where pairs of
    sites with r^2 <= threshold are True
    """
    print("Calculating LD arrays... (Warning: Very slow!)")
    mask=np.zeros(ts.num_sites,dtype="bool")

    def thread_worker(thread_index):
        ld = msp.LdCalculator(ts)
        chunk_size = int(math.ceil(len(mask) / num_threads))
        nextSite = thread_index * chunk_size
        stop = nextSite + chunk_size
        while True:
            mask[nextSite] = True
            r2 = (ld.r2_array(nextSite) <= thresh)
            if nextSite > stop or len(r2) == 0 or not np.any(r2):
                break
            nextSite += (1 + np.argmax(r2))

    threads = [
            threading.Thread(target=thread_worker, args=(j,))
            for j in range(num_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print("LD calculation finished!")
    return mask


def simulate(out_path, species, model, genetic_map, seed, chrmStr,
             sample_size=20, population=0, ld_thresh=1.0, max_workers=1):
    mask_path = out_path+".r2Mask.p"
    sfs_path = out_path+".sfs.pdf"
    chrom = species.genome.chromosomes[chrmStr]
    samples = [msp.Sample(population=population, time=0)] * sample_size
    print("Simulating...")
    ts = msp.simulate(
        samples=samples,
        recombination_map=chrom.recombination_map(genetic_map.name),
        mutation_rate=chrom.default_mutation_rate,
        random_seed=seed,
        **model.asdict())
    ts.dump(out_path)
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    SFSs=[]
    SFSs.append(allel.sfs(haps.count_alleles()[:,1])[1:])
    print("Simulation finished!")
    if ld_thresh < 1.0:
        ul = unlinked(ts,ld_thresh,max_workers)
        mask_file = open(mask_path, "wb")
        pickle.dump(ul, mask_file)
        SFSs.append(allel.sfs(haps[ul,:].count_alleles()[:,1])[1:])
    plot_sfs(SFSs,sfs_path)
