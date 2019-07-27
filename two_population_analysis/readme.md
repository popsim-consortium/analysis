# Readme for two population analysis pipeline

This directory contains the code to run analysis of
demographic inference using dadi and fastsimcoal2 programs, with
identical resulting data from any demographic model
found in `stdpopsim` as input to each.
This provides an example of how standardized
population simulations may be used to compare methods.

The Snakemake workflow includes the necessary pipeline
for simulation, analyses, and plotting in an efficient manner.
Simply choose your parameters in a config file,
and let snakemake handle the rest.
(For large runs, the use of a cluster is highly encouraged)

_NOTE_: This requires a bleeding edge install of `msprime` at the
moment.


## Dependencies
Dadi needs to be installed before running this pipeline: https://bitbucket.org/gutenkunstlab/dadi/wiki/Installation  

To install dadi, first install scipy and numpy to the popsim_env_test conda environment:
```
conda install numpy
conda install scipy
```
Now we can install dadi:
```
git clone https://bitbucket.org/gutenkunstlab/dadi
cd dadi
python setup.py install
cd ..
```
To install fastsimcoal2, download compiled version for your platform here: http://cmpg.unibe.ch/software/fastsimcoal2  
Make sure to download into the two_population directory.
Next, we unzip the download and rename the fastsimcoal directory:  
```
unzip fsc26_*.zip
rm fsc26_*.zip
mv fsc26_* fsc26
```
For mac and linux users, you will also need to change the permissions on the fsc code:
```
cd fsc26/
chmod +x fsc26
cd ..
```
And we can install smc++ using this command
```
conda install -c terhorst -c bioconda smcpp
```
Finally, mac users may also have to update their gcc if you do not have a recent version installed:
```
wget http://prdownloads.sourceforge.net/hpc/gcc-7.1-bin.tar.gz
gunzip gcc-7.1-bin.tar.gz
sudo tar -xvf gcc-7.1-bin.tar -C /.
```


## Workflow

The analysis includes two programs for inferring population size, split times,
and migration rates from two populations:
[dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master/), and
[fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/).

To run an analysis, create a directory (wherever you want)
where all results, and intermediate
files will be stored. Next, create and place a file named `config.json` in it.
The json file must contain key : value combos described below. An example
might look like this:

```json
{
    "seed" : 12345,
    "num_sampled_genomes_per_replicate" : 20,
    "replicates" : 10,
    "species" : "homo_sapiens",
    "model" : "GutenkunstThreePopOutOfAfrica",
    "genetic_map" : "HapmapII_GRCh37",
    "chrm_list" : "chr21,chr22",
}
```

Once you have creates a directory which contains the config file
simply run snakemake from _this_ directory (two_population_analysis), and point it to your analysis run
directory, like so

`snakemake -j 40 --config config="/projects/kernlab/jgallowa/homo_sapiens_Gutenkunst_0"`

where `-j` is the number cores available to run jobs in parallel, and
`--config` points to the _directory_ that contains the config file.


### Cluster environments
Our workflow can also be run on a cluster. To do so requires
the setup of a `.json` configuration file that lets `snakemake`
know about your cluster. We have provided an example of
such a file in `cluster_talapas.json` that is for use with a
University of Oregon SLURM cluster.

```json
{
    "__default__" :
    {
        "time" : "02:30:00",
        "n" : 1,
        "partition" : "kern,preempt",
        "mem" : "16G",
        "cores" : "4",
    },
    "run_msmc" :
    {
        "time" : "05:00:00",
        "mem" : "32G",
        "cores" : "4",
    }

}
```

At a minimum, you should
edit the names of the partition to match those on your own HPC.
The workflow can then be launched with the call

`snakemake -j 999 --config config="/projects/kernlab/jgallowa/homo_sapiens_Gutenkunst_0" --cluster-config cluster_talapas.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} --mem-per-cpu {cluster.mem} -c {cluster.cores}"`

(it may prove useful to simply put this command in a bash script)

and jobs will be automatically farmed out to the cluster

### Final output
The current final output is are two plots: one comparing population size  
and divergence time estimates, and a second comparing migration rate estimates, i.e.,  
`homo_sapiens_Gutenkunst/estimates_mig_dadi_fsc.png`  
`homo_sapiens_Gutenkunst/estimates_N_tdiv_dadi_fsc.png`


## Parameter Description

`seed` : `<class 'int'>`
This sets the seed such that any anaysis configuration
and run can be replicated exactly.

`num_sampled_genomes_per_replicate` : `<class 'int'>`
This is the haploid number
of genomes to simulate for an analysis run.

`replicates` : `<class 'int'>` The number of replicate simulations to run and
analyze.

`output_dir` : `<class 'str'>` This where all the intermediate files and results
will be stored. This includes all plots preduced, as well as all input/output files
from the simulations and software runs organized into subdirectories by
replicate seed.

`species` : `<class 'module'>` from `stdpopsim` such as `stdpopsim.homo_sapiens`.
Chromosomes, demographics, and recombination maps should be derived from this.

`model` : `<class childclass Model>`
This is the class that defines the demography associated with a species. All models
should inherit from `<class Model>`

`genetic_map` : `<class childclass GeneticMap>` This will define the genetic map
used for simulations.

`chrm_list` : `<class 'str'>` A string of the chromosome names you would like to simulate,
separated by commas. All chromosomes simulated will be fed
as a single input into each analysis by the inference programs, for each replicate.
Set to "all" to simulate all chromsomes for the genome.
