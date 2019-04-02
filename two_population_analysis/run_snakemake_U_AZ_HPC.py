#!/usr/bin/env python

# Name of job
#PBS -N run_snakemake

# Maximum CPU and clock time per job
#PBS -l cput=1:00:0
#PBS -l walltime=1:00:0

# Submit an array of jobs
###PBS -J 1-100
# Which queue should we use?
#PBS -q oc_high_pri
###PBS -q oc_standard
###PBS -q oc_windfall

# Ensure your time gets charged to the group
#PBS -W group_list=rgutenk
# Pass your environment variables to the pbs job
#PBS -V
# Join standard out and standar error messages
#PBS -j oe


#
# For a multiprocessor job that needs to grab a whole node
#
#PBS -l select=1:ncpus=1:mem=6gb:pcmem=6gb
#PBS -l place=pack:shared

import sys, os

workdir = os.environ.get('PBS_O_WORKDIR')
print(workdir)
try:
	os.chdir(workdir)
except:
	pass

os.system('snakemake -j 200 --cluster-config U_AZ_hpc.json --cluster "qsub -N {cluster.N} -W {cluster.W} -q {cluster.q} -l {cluster.cpu_setup} -l {cluster.cput} -l {cluster.walltime} -j oe -V"')


