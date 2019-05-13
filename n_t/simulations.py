"""
Module to run the required simulations.
"""

import msprime
from stdpopsim import homo_sapiens


def simulate(out_path, species, model, genetic_map, seed, chrmStr, 
             sample_size=20, population=0):
    chrom = species.genome.chromosomes[chrmStr]
    # TODO : Sample outside of population 0?
    samples = [msprime.Sample(population=population, time=0)] * sample_size
    ts = msprime.simulate(
        samples=samples,
        recombination_map=chrom.recombination_map(genetic_map.name),
        mutation_rate=chrom.default_mutation_rate,
        random_seed=seed,
        **model.asdict())
    ts.dump(out_path)
