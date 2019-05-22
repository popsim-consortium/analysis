"""
Utilities for working with msmc
"""
import logging
import subprocess
import tskit
import numpy as np
import os

def prune_tree_sequence(tree_sequence_path, num_samples):
    """
    take in a tree sequence, and a number of samples
    less than the number of samples in the tree, 
    then simplify the tree sequence on that subset.
    """
    ts = tskit.load(tree_sequence_path)
    if ts.num_samples > num_samples:
        subset = np.random.choice(ts.samples(), num_samples, replace=False)
        ts = ts.simplify(subset)
    return ts


def write_msmc_file(path, num_sampled_genomes_msmc):
    """
    take one .trees file and write out 
    path.multihep.txt which acts as a single input to msmc

    This seems hacky atm, but let's getting working then hash out
    the details
    """
    for sample_size in num_sampled_genomes_msmc:
        ts = prune_tree_sequence(path, sample_size)
        dirr = os.path.dirname(path)
        filen = os.path.basename(path)
        sep = filen.split(".")
        chrom = sep[0]
        sep.insert(0,str(sample_size))
        output = os.path.join(dirr,".".join(sep) + ".multihep.txt")
        fi = open(output, "w")
        prev = 0
        for var in ts.variants():
            cur = int(var.site.position)
            if cur > prev:
                geno = ''.join(map(str,var.genotypes))
                fi.write(f"{chrom}\t{cur}\t{cur-prev}\t{geno}\n")
            prev = cur
        fi.close()
    return None



def run_msmc_estimate(input_files, output_file, msmc_exec_loc, iterations=1, ncores=1):
    """
    This is to run the msmc command and get estimates,
    It then will convert the scales times and pop sizes
    into generation years and population sizes.

    The final estimates will be written to 
    input_file.final.txt
    """
    
    cmd = (f"{msmc_exec_loc} --fixedRecombination -o \
        {output_file} -i {iterations} {input_files}")
    subprocess.run(cmd, shell=True, check=True)
    return None


def convert_msmc_output(results_file, outfile, mutation_rate, generation_time):
    """
    This function converts the output from msmc into a csv the will be read in for
    plotting comparison.
    
    MSMC outputs times and rates scaled by the mutation rate per basepair per generation. 
    First, scaled times are given in units of the per-generation mutation rate. 
    This means that in order to convert scaled times to generations, 
    divide them by the mutation rate. In humans, we used mu=1e-8 per basepair per generation.
    To convert generations into years, multiply by the generation time, for which we used 10 years.

    To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, 
    scaledPopSize = 1 / lambda00. Then divide this scaled population size by 2*mu
    """
    out_fp = open(outfile, "w")
    in_fp = open(results_file, "r")
    header = in_fp.readline()
    out_fp.write("label,x,y,plot_type,plot_num")
    for line in in_fp:
        result = line.split()
        time = float(result[1])
        time_generation = time / mutation_rate
        time_years = time_generation * generation_time
        lambda00 = float(result[3])
        scaled_pop_size = 1 / lambda00
        size = scaled_pop_size / (2*mutation_rate)
        out_fp.write(f"pop0,{time_years},{size},path,0\n")
    out_fp.close
    return None
