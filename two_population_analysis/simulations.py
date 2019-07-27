"""
Module to run the required simulations.
"""

import msprime
from stdpopsim import homo_sapiens


def simulate(out_path, species, model, genetic_map, seed, chrmStr, sample_size=20):
    chrom = species.genome.chromosomes[chrmStr]
    samples_pops_joint = [msprime.Sample(population=0, time=0)] * sample_size + [msprime.Sample(population=1, time=0)] * sample_size
    ts_pops_joint = msprime.simulate(
    	samples=samples_pops_joint,
    	recombination_map=chrom.recombination_map(),
    	mutation_rate=chrom.default_mutation_rate,
    	random_seed=seed,
    	**model.asdict())
    ts_pops_joint.dump(out_path)
