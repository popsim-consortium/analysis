# analysis
Analysis of inference methods on standard population models

# how to set up your environment to run the analysis
We recommend you start by creating a new `conda` environment for the analysis. 

```
conda create -n popsim_env_test python=3.6.8 --yes
conda activate popsim_env_test
```

Next clone and install `stdpopsim`
```
git clone https://github.com/popgensims/stdpopsim.git
cd stdpopsim
python setup.py install
```

Now clone the analysis repo, and install its dependencies
```
git clone https://github.com/popgensims/analysis.git
cd analysis/

for c in terhorst bioconda defaults conda-forge; do conda config --add channels $c; done
conda install --file requirements.txt --yes
````

Further instructions can be currently found in each task directory
