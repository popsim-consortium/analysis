#Readme for N(t) example
This directory has a quick example of a workflow for estimating changing 
population size over time from a simulated sample ($N(t)$) using
stairwayplot. The basic workflow is organized using `snakemake`
so the only step needed to perform the simulation and run 
the estimation is

`$ snakemake -j 8`

In this case the `-j 8` option tells `snakemake` to use 8 cores for 
the workflow. If more cores are available by all means increase this 
number.

If something goes wrong or you need to start from scratch I've
included a `clean` rule in the `Snakefile`, so you can clean 
everything out using

`$ snakemake clean`

###Cluster environments
Our workflow can also be run on a cluster. To do so requires
the setup of a `.json` configuration file that lets `snakemake`
know about your cluster. We have provided an example of 
such a file in `cluster_talapas.json` that is for use with a
University of Oregon SLURM cluster. At a minimum, you should
edit the names of the partition to match those on your own HPC.
The workflow can then be launched with the call

`$ snakemake -j 999 --cluster-config cluster_talapas.json --cluster "sbatch -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"

and jobs will be automatically farmed out to the cluster
