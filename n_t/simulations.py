"""
Module to run the required simulations
and create unlinked-sites mask
"""
import msprime as msp


def simulate(out_path, species, model, genetic_map, seed, chrmStr,
             sample_size=20, population=0):
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
    print("Simulation finished!")
