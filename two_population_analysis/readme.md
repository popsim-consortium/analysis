# Readme for two population analysis

This directory contains a snakemake pipeline for running dadi and fastsimcoal2 on data simulated with msprime.

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
Finally, mac users may also have to update their gcc if you do not have a recent version installed:
```
wget http://prdownloads.sourceforge.net/hpc/gcc-7.1-bin.tar.gz
gunzip gcc-7.1-bin.tar.gz
sudo tar -xvf gcc-7.1-bin.tar -C /.
```
