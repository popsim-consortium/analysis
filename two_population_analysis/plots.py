"""
Code for generating plots.
"""
import pandas
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

sns.set_style("darkgrid")


def plot_fsc_dadi_results(dadi_infile,fsc_infile, outfile, simulated_genome_length):


    ## convert dadi and fsc estimates to be comparable and add "A-" needed for converting df to long format
    dadi = pandas.read_csv(dadi_infile, skiprows=0, sep="\t")
    dadi['A-Na'] = dadi['theta']/(4*1e-8*simulated_genome_length)
    dadi['A-N2'] = dadi['nu1']*dadi['A-Na']
    dadi['A-N1'] = dadi['nu2']*dadi['A-Na']
    dadi['A-TDIV'] = dadi['T']*dadi['A-Na']*2
    dadi['A-MIG12'] = dadi['m1']/(dadi['A-Na']*2)
    dadi['A-MIG21'] = dadi['m2']/(dadi['A-Na']*2)

    fsc = pandas.read_csv(fsc_infile, skiprows=0, sep="\t")
    fsc['A-Na'] = fsc['ANCSIZE']/2
    fsc['A-N1'] = fsc['NPOP1']/2
    fsc['A-N2'] = fsc['NPOP2']/2
    fsc['A-TDIV'] = fsc['TDIV']
    fsc['A-MIG12'] =fsc['MIG12']
    fsc['A-MIG21'] =fsc['MIG21']


    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(7, 7))
    truth_N = pandas.DataFrame(np.array([[7300,12300,3517.113,140000/25]]), columns=['A-Na','A-N2','A-N1','A-TDIV'])

    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_N=data_N.reset_index(drop=True).reset_index()
    data_N_long=pandas.wide_to_long(df=data_N,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_N_long,x='parameter',y='estimate', hue='method', jitter=False, palette='muted') 
    f.savefig(outfile[0], bbox_inches='tight', alpha=0.8)


    #create and plot data frame of migration rate estimates

    f, ax = plt.subplots(figsize=(7, 7))
    truth_M = pandas.DataFrame(np.array([[3e-5,3e-5]]), columns=['A-MIG12','A-MIG21'])

    data_M = fsc[['A-MIG12','A-MIG21']]
    data_M=data_M.append(dadi[['A-MIG12','A-MIG21']], sort=True)
    data_M=data_M.append(truth_M, sort=True)
    num_sims = dadi.shape[0]
    data_M['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_M=data_M.reset_index(drop=True).reset_index()
    data_M_long=pandas.wide_to_long(df=data_M,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_M_long,x='parameter',y='estimate', hue='method', jitter=False, palette='muted')
    ax.set(ylim=(0, 4e-4))
    f.savefig(outfile[1], bbox_inches='tight', alpha=0.8)
