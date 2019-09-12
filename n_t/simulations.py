"""
Module to run the required simulations
and create unlinked-sites mask
"""
import msprime as msp


def simulate(out_path, species, model, genetic_map, seed, chrmStr,
             sample_size=20, population=0):
    chrom = species.genome.chromosomes[chrmStr]
    samples = [msp.Sample(population=population, time=0)] * sample_size
    print(chrmStr)
    print("Simulating...")
    if genetic_map:
        ts = msp.simulate(
            samples=samples,
            recombination_map=chrom.recombination_map(genetic_map.name),
            mutation_rate=chrom.default_mutation_rate,
            random_seed=seed,
            **model.asdict())
    else:
        print(f"log check: {chrom.default_recombination_rate}")
        ts = msp.simulate(
            samples=samples,
            recombination_map=msp.RecombinationMap.uniform_map(
                chrom.length, chrom.default_recombination_rate),
            mutation_rate=chrom.default_mutation_rate,
            random_seed=seed,
            **model.asdict())
    ts.dump(out_path)
    print("Simulation finished!")
