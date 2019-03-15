"""
Module to run the required simulations.
"""

import msprime
from stdpopsim import homo_sapiens


def homo_sapiens_Gutenkunst(path, seed, chrmStr, sample_size=20):
    chrom = homo_sapiens.genome.chromosomes[chrmStr]

    model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
    # model.debug()

    # Currently sampling 20 individuals from a single popn.
    samples = [msprime.Sample(population=0, time=0)] * sample_size
    ts = msprime.simulate(
        samples=samples,
        recombination_map=chrom.recombination_map(),
        mutation_rate=chrom.mean_mutation_rate,
        random_seed=seed,
        **model.asdict())
    ts.dump(path)
