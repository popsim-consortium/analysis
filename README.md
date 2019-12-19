# analysis
Analysis of inference methods on standard population models

# how to set up your environment to run the analysis
We recommend you start by creating a new `conda` environment for the analysis. 

```
conda create -n popsim_env_test --yes -c conda-forge -c terhorst smcpp python=3.7 
conda activate popsim_env_test
```

Next, install `stdpopsim`
```
python3 -m pip install stdpopsim==0.1.0
```

Now clone the analysis repo, and install its dependencies
```
git clone https://github.com/popgensims/analysis.git
cd analysis/

conda install --file requirements.txt --yes
````

For using `msmc` we need to download and compile it to play nice
with the conda environment that we have set up.
```
cd extern
git clone https://github.com/stschiff/msmc.git
cat msmc_makefile_stdpopsim_patch > msmc/Makefile && cd msmc && make
cd ../../
```

Further instructions can be currently found in each task directory
