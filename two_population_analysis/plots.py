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


def plot_fsc_dadi_smcpp_results_human_IM(dadi_infile,fsc_infile,smcpp_infile, outfile, simulated_genome_length):

    sns.set(style="whitegrid")

    ## convert dadi and fsc estimates to be comparable and add "A-" needed for converting df to long format
    dadi = pandas.read_csv(dadi_infile, skiprows=0, sep="\t")
    dadi['A-Na'] = dadi['theta']/(4*1e-8*simulated_genome_length)
    dadi['A-N2'] = dadi['nu2']*dadi['A-Na']
    dadi['A-N1'] = dadi['nu1']*dadi['A-Na']
    dadi['A-TDIV'] = dadi['T']*dadi['A-Na']*2
    dadi['A-MIG12'] = dadi['m1']/(dadi['A-Na']*2)
    dadi['A-MIG21'] = dadi['m2']/(dadi['A-Na']*2)

    fsc = pandas.read_csv(fsc_infile, skiprows=0, sep="\t")
    fsc['A-Na'] = fsc['ANCSIZE']/2
    fsc['A-N1'] = fsc['NPOP2']/2
    fsc['A-N2'] = fsc['NPOP1']/2
    fsc['A-TDIV'] = fsc['TDIV']
    fsc['A-MIG12'] =fsc['MIG12']
    fsc['A-MIG21'] =fsc['MIG21']

    smcpp = pandas.read_csv(smcpp_infile, skiprows=0, sep="\t")
    smcpp['A-N1'] = smcpp['HarmonicMeanPop1']
    smcpp['A-N2'] = smcpp['HarmonicMeanPop2']
    smcpp['A-Na'] = None
    smcpp['A-TDIV'] = smcpp['DivTime']


    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(10, 10))
    # need to convert Tdiv from years to generations
    # Neu (N1) computed as harmonic mean over history of european pop
    truth_N = pandas.DataFrame(np.array([[7300,12300,3517.113,140000/25]]), columns=['A-Na','A-N1','A-N2','A-TDIV'])

    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(smcpp[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['smcpp']*num_sims+['truth']
    data_N=data_N.reset_index(drop=True).reset_index()
    data_N_long=pandas.wide_to_long(df=data_N,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    # add 'truth' from mean coalescence time method
    meal_coal_time_N1 = pandas.DataFrame(np.array([["N1","N_mean_coal_time",18114.993411358824/2]]), columns=['parameter','method','estimate'])
    meal_coal_time_N2 = pandas.DataFrame(np.array([["N2","N_mean_coal_time",13021.533578633795/2]]), columns=['parameter','method','estimate'])
    data_N_long = data_N_long.append(meal_coal_time_N1,sort=True)
    data_N_long = data_N_long.append(meal_coal_time_N2,sort=True)

    # add 'truth' from: N = - T/(2*log(P(not coalesced at split)))
    truth_p_not_coal_N1 = pandas.DataFrame(np.array([["N1","N_p_not_coal",11671.97129362]]), columns=['parameter','method','estimate'])
    truth_p_not_coal_N2 = pandas.DataFrame(np.array([["N2","N_p_not_coal",5352.07636014]]), columns=['parameter','method','estimate'])
    data_N_long = data_N_long.append(truth_p_not_coal_N1,sort=True)
    data_N_long = data_N_long.append(truth_p_not_coal_N2,sort=True)


    sns.stripplot(data=data_N_long,x='parameter',y='estimate', hue='method', jitter=True, palette='muted', size=7, linewidth=1).set_title('Human IM inference')
    ax.set(ylim=(0, 20000))
    f.savefig(outfile[0], bbox_inches='tight', alpha=0.8)


    #create and plot data frame of migration rate estimates

    f, ax = plt.subplots(figsize=(10, 10))
    truth_M = pandas.DataFrame(np.array([[3e-5,3e-5]]), columns=['A-MIG12','A-MIG21'])

    data_M = fsc[['A-MIG12','A-MIG21']]
    data_M=data_M.append(dadi[['A-MIG12','A-MIG21']], sort=True)
    data_M=data_M.append(truth_M, sort=True)
    num_sims = dadi.shape[0]
    data_M['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_M=data_M.reset_index(drop=True).reset_index()
    data_M_long=pandas.wide_to_long(df=data_M,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_M_long,x='parameter',y='estimate', hue='method', jitter=True, palette='muted', size=7, linewidth=1)
    ax.set(ylim=(-1e-5, 4e-4))
    f.savefig(outfile[1], bbox_inches='tight', alpha=0.8)




def plot_fsc_dadi_smcpp_results_drosophila_IM(dadi_infile,fsc_infile, smcpp_infile, outfile, simulated_genome_length):

    ## convert dadi and fsc estimates to be comparable and add "A-" needed for converting df to long format
    dadi = pandas.read_csv(dadi_infile, skiprows=0, sep="\t")
    dadi['A-Na'] = dadi['theta']/(4*8.4e-9*simulated_genome_length)
    dadi['A-N2'] = dadi['nu2']*dadi['A-Na']
    dadi['A-N1'] = dadi['nu1']*dadi['A-Na']
    dadi['A-TDIV'] = dadi['T']*dadi['A-Na']*2
    dadi['A-MIG12'] = dadi['m1']/(dadi['A-Na']*2)
    dadi['A-MIG21'] = dadi['m2']/(dadi['A-Na']*2)

    fsc = pandas.read_csv(fsc_infile, skiprows=0, sep="\t")
    fsc['A-Na'] = fsc['ANCSIZE']/2
    fsc['A-N1'] = fsc['NPOP2']/2
    fsc['A-N2'] = fsc['NPOP1']/2
    fsc['A-TDIV'] = fsc['TDIV']
    fsc['A-MIG12'] =fsc['MIG12']
    fsc['A-MIG21'] =fsc['MIG21']

    smcpp = pandas.read_csv(smcpp_infile, skiprows=0, sep="\t")
    smcpp['A-N1'] = smcpp['HarmonicMeanPop1']
    smcpp['A-N2'] = smcpp['HarmonicMeanPop2']
    smcpp['A-Na'] = None
    smcpp['A-TDIV'] = smcpp['DivTime']


    sns.set(style="whitegrid")
    #create and plot data frame of pop size and t div estimates
    f, ax = plt.subplots(figsize=(10, 10))

    # value for Neu (N1) calculated as harmonic mean of Ne1 and Ne0
    truth_N = pandas.DataFrame(np.array([[8.603e06,8.603e06,95365.76,158000]]), columns=['A-Na','A-N1','A-N2','A-TDIV'])
    #truth_N1 = pandas.DataFrame(np.array([[8.603e06/5,-1,-1,-1]]), columns=['A-Na','A-N2','A-N1','A-TDIV'])


    data_N = fsc[['A-Na','A-N1','A-N2','A-TDIV']]
    data_N=data_N.append(dadi[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(smcpp[['A-Na','A-N1','A-N2','A-TDIV']], sort=True)
    data_N=data_N.append(truth_N, sort=True)
    #data_N=data_N.append(truth_N1, sort=True)

    num_sims = dadi.shape[0]
    data_N['method']=['fsc']*num_sims+['dadi']*num_sims+['smcpp']*num_sims+['truth']
    data_N=data_N.reset_index(drop=True).reset_index()
    data_N_long=pandas.wide_to_long(df=data_N,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    #add extra Na point for Li & Stephan model
    Na_truth = pandas.DataFrame(np.array([["Na","truth",8.603e06/5]]), columns=['parameter','method','estimate'])
    data_N_long = data_N_long.append(Na_truth,sort=True)

    # add 'truth' from mean coalescence time method
    mean_coal_time_N1 = pandas.DataFrame(np.array([["N1","N_mean_coal_time",3910531.1902434235/2]]), columns=['parameter','method','estimate'])
    mean_coal_time_N2 = pandas.DataFrame(np.array([["N2","N_mean_coal_time",1779082.8243646843/2]]), columns=['parameter','method','estimate'])
    data_N_long = data_N_long.append(mean_coal_time_N1,sort=True)
    data_N_long = data_N_long.append(mean_coal_time_N2,sort=True)

    # add 'truth' from: N = - T/(2*log(P(not coalesced at split)))
    truth_p_not_coal_N1 = pandas.DataFrame(np.array([["N1","N_p_not_coal",8603000.00000001]]), columns=['parameter','method','estimate'])
    truth_p_not_coal_N2 = pandas.DataFrame(np.array([["N2","N_p_not_coal",344699.62629219]]), columns=['parameter','method','estimate'])
    data_N_long = data_N_long.append(truth_p_not_coal_N1,sort=True)
    data_N_long = data_N_long.append(truth_p_not_coal_N2,sort=True)


    sns.stripplot(data=data_N_long,x='parameter',y='estimate', hue='method', jitter=True, palette='muted', size=7, linewidth=1)
    f.savefig(outfile[0], bbox_inches='tight', alpha=0.8)


    #create and plot data frame of migration rate estimates

    f, ax = plt.subplots(figsize=(10, 10))
    truth_M = pandas.DataFrame(np.array([[0,0]]), columns=['A-MIG12','A-MIG21'])

    data_M = fsc[['A-MIG12','A-MIG21']]
    data_M=data_M.append(dadi[['A-MIG12','A-MIG21']], sort=True)
    data_M=data_M.append(truth_M, sort=True)
    num_sims = dadi.shape[0]
    data_M['method']=['fsc']*num_sims+['dadi']*num_sims+['truth']
    data_M=data_M.reset_index(drop=True).reset_index()
    data_M_long=pandas.wide_to_long(df=data_M,stubnames=["A"],i='index',j='parameter',sep="-",suffix='(\d+|\w+)').reset_index().rename(columns={'A':'estimate','method':'method'})

    sns.stripplot(data=data_M_long,x='parameter',y='estimate', hue='method', jitter=True, palette='muted', size=7, linewidth=1)
    ax.set(ylim=(-1e-7, 1e-6))
    f.savefig(outfile[1], bbox_inches='tight', alpha=0.8)





def plot_compound_smcpp(infiles, outfile):
    cmd = f"smc++ plot {outfile} {infiles}"
    #for infile in infiles:
    #cmd = cmd = + f" {infile}"
    #logging.info("Running:" + cmd)
    print ("Running this! " + cmd)
    subprocess.run(cmd, shell=True, check=True)
