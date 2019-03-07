"""
Example simulation for analysis
"""
import msprime


# import stdpopsim
from stdpopsim import homo_sapiens

chrom = homo_sapiens.genome.chromosomes["chr22"]
recomb_map = chrom.recombination_map()

model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
model.debug()

# Currently sampling 20 individuals from a single popn.
tmp_samples = [
    msprime.Sample(population=0, time=0)
]
samples = tmp_samples * 20
ts = msprime.simulate(
    samples=samples,
    recombination_map=chrom.recombination_map(),
    mutation_rate=chrom.mean_mutation_rate,
    **model.asdict())
# Hard coded output name. FIX ME
ts.dump("G3OOA.trees")
