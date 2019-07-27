"""
Code for generating plots.
"""
import pandas
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import subprocess

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

sns.set_style("darkgrid")


def plot_fsc_dadi_results_human_IM(dadi_infile,fsc_infile, outfile, simulated_genome_length):
    
    
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





def plot_fsc_dadi_results_drosophila_IM(dadi_infile,fsc_infile, outfile, simulated_genome_length):
    
    
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
    
    
    smcpp = pandas.read_csv(smcpp_infile, skiprows=0, sep="\t")
    smcpp['A-N1'] = smcpp['HarmonicMeanPop1']
    smcpp['A-N2'] = smcpp['HarmonicMeanPop2']
    smcpp['A-Na'] = None
    smcpp['A-TDIV'] = smcpp['DivTime']
    
    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(7, 7))
    
    truth_N = pandas.DataFrame(np.array([[5e6,8.603e06,1051914,1580]]), columns=['A-Na','A-N2','A-N1','A-TDIV'])
    
    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(smcpp[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['smcpp']*num_sims+['truth']
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


def plot_fsc_dadi_smcpp_results_human_IM(dadi_infile,fsc_infile, smcpp_infile, outfile, simulated_genome_length):
    
    
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
    
    
    smcpp = pandas.read_csv(smcpp_infile, skiprows=0, sep="\t")
    smcpp['A-N1'] = smcpp['HarmonicMeanPop1']
    smcpp['A-N2'] = smcpp['HarmonicMeanPop2']
    smcpp['A-Na'] = None
    smcpp['A-TDIV'] = smcpp['DivTime']
    
    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(7, 7))
    truth_N = pandas.DataFrame(np.array([[7300,12300,3517.113,140000/25]]), columns=['A-Na','A-N2','A-N1','A-TDIV'])
    
    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(smcpp[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['smcpp']*num_sims+['truth']
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





def plot_fsc_dadi_smcpp_results_drosophila_IM(dadi_infile,fsc_infile, outfile, simulated_genome_length):
    
    ## convert dadi and fsc estimates to be comparable and add "A-" needed for converting df to long format
    dadi = pandas.read_csv(dadi_infile, skiprows=0, sep="\t")
    dadi['A-Na'] = dadi['theta']/(4*8.4e-9*simulated_genome_length)
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

    smcpp = pandas.read_csv(fsc_infile, skiprows=0, sep="\t")
    smcpp['A-N1'] = fsc['HarmonicMeanPop1']
    smcpp['A-N2'] = fsc['HarmonicMeanPop2']
    smcpp['A-Na'] = None
    smcpp['A-TDIV'] = fsc['DivTime']


    #create and plot data frame of pop size and t div estimates
    # TODO: what are the ground truth parameter values?
    f, ax = plt.subplots(figsize=(7, 7))
    truth_N = pandas.DataFrame(np.array([[5e6,8.603e06,1051914,1580]]), columns=['A-Na','A-N2','A-N1','A-TDIV'])
    
    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(smcpp[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['smcpp']*num_sims+['truth']
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







def plot_fsc_dadi_results_gute2pop(dadi_infile,fsc_infile, outfile, simulated_genome_length):


    ## convert dadi and fsc estimates to be comparable and add "A-" needed for converting df to long format
    dadi = pandas.read_csv(dadi_infile, skiprows=0, sep="\t")
    dadi['A-NANC'] = dadi['theta']/(4*1e-8*simulated_genome_length)
    dadi['A-NEU'] = dadi['nuEu']*dadi['A-NANC']
    dadi['A-NAF'] = dadi['nuAf']*dadi['A-NANC']
    dadi['A-NBOT'] = dadi['nuB']*dadi['A-NANC']
    dadi['A-TANC'] = dadi['TAf']*dadi['A-NANC']*2
    dadi['A-TEU'] = dadi['TEuF']*dadi['A-NANC']*2
    dadi['A-TOOA'] = dadi['TB']*dadi['A-NANC']*2
    dadi['A-MIG_AFEU'] = dadi['mAfEu']/(dadi['A-NANC']*2)
    dadi['A-MIG_AFB'] = dadi['mAfB']/(dadi['A-NANC']*2)

    fsc = pandas.read_csv(fsc_infile, skiprows=0, sep="\t")
    fsc['A-NANC'] = fsc['NANC']/2
    fsc['A-NEU'] = fsc['NEU']/2
    fsc['A-NAF'] = fsc['NAF']/2
    fsc['A-NBOT'] = fsc['NBOT']/2
    fsc['A-TEU'] =fsc['TEU']
    fsc['A-TOOA'] =fsc['TOOA']
    fsc['A-TANC'] =fsc['TANC']
    fsc['A-MIG_AFEU'] =fsc['MIG_AFEU']
    fsc['A-MIG_AFB'] =fsc['MIG_AFB']


    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(7, 7))
    truth_N = pandas.DataFrame(np.array([[7300,29725.34,12300,2100,21.2e3/25,140e3/25,220e3/25]]), columns=['A-NANC','A-NEU','A-NAF','A-NBOT','A-TEU','A-TOOA','A-TANC'])

    data_N = fsc[['A-NANC','A-NEU','A-NAF','A-NBOT','A-TEU','A-TOOA','A-TANC']]
    data_N=data_N.append(dadi[['A-NANC','A-NEU','A-NAF','A-NBOT','A-TEU','A-TOOA','A-TANC']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_N=data_N.reset_index(drop=True).reset_index()
    data_N_long=pandas.wide_to_long(df=data_N,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_N_long,x='parameter',y='estimate', hue='method', jitter=False, palette='muted')
    f.savefig(outfile[0], bbox_inches='tight', alpha=0.8)


    #create and plot data frame of migration rate estimates

    f, ax = plt.subplots(figsize=(7, 7))
    truth_M = pandas.DataFrame(np.array([[25e-5,3e-5]]), columns=['A-MIG_AFB','A-MIG_AFEU'])

    data_M = fsc[['A-MIG_AFB','A-MIG_AFEU']]
    data_M=data_M.append(dadi[['A-MIG_AFB','A-MIG_AFEU']], sort=True)
    data_M=data_M.append(truth_M, sort=True)
    num_sims = dadi.shape[0]
    data_M['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_M=data_M.reset_index(drop=True).reset_index()
    data_M_long=pandas.wide_to_long(df=data_M,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_M_long,x='parameter',y='estimate', hue='method', jitter=False, palette='muted')
    ax.set(ylim=(-1e-4, 4e-4))
    f.savefig(outfile[1], bbox_inches='tight', alpha=0.8)




def plot_compound_smcpp(infiles, outfile):
    cmd = f"smc++ plot {outfile} {infiles}"
    #for infile in infiles:
    #cmd = cmd = + f" {infile}"
    #logging.info("Running:" + cmd)
    print ("Running this! " + cmd)
    subprocess.run(cmd, shell=True, check=True)



