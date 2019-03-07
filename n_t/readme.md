#Readme for N(t) example
This directory has a quick example of a workflow for estimating changing 
population size over time from a simulated sample ($N(t)$) using
stairwayplot. The basic workflow is organized using `snakemake`
so the only step needed to perform the simulation and run 
the estimation is

`$ snakemake`

If something goes wrong or you need to start from scratch I've
included a `clean` rule in the `Snakefile`, so you can clean 
everything out using

`$ snakemake clean`
