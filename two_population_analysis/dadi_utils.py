
import allel
import numpy as np
import stdpopsim
import dadi
import os
import matplotlib.pyplot as plt
import random
import tskit
import pandas as pd




def ts_to_dadi_sfs(ts_path,out_path,out_path_nonvariant, sample_size=20, mask_file=None):
	'''
	Generate however many different SFS with msprime and convert+save them into SFS for dadi to use.
	'''
	ts = tskit.load(ts_path)

	#haps_pops_joint = np.array(ts.genotype_matrix())

	haps = ts.genotype_matrix()

	total_length = ts.sequence_length

	# Masking
	retain = np.full(ts.get_num_mutations(), False)
	if mask_file:
		mask_table = pd.read_csv(mask_file, sep="\t", header=None)
		chrom = ts_path.split("/")[-1].split(".")[0]
		sub = mask_table[mask_table[0] == chrom]
		mask_ints = pd.IntervalIndex.from_arrays(sub[1], sub[2])
		snp_locs = [int(x.site.position) for x in ts.variants()]
		tmp_bool = [mask_ints.contains(x) for x in snp_locs]
		retain = np.logical_or(retain, tmp_bool)
		#print(retain)
		total_length -= np.sum(mask_ints.length)
	#print(ts.sequence_length)
	#print(total_length)

	retain = np.logical_not(retain)

	haps_pops_joint= np.array(haps[retain, :])

	#Break up the haplotypes into seperate populations based on sample_size
	haps_pop0_joint = haps_pops_joint[:,:sample_size]
	haps_pop1_joint = haps_pops_joint[:,sample_size:]


	genotypes_pop0_joint = allel.HaplotypeArray(haps_pop0_joint).to_genotypes(ploidy=2)
	allele_counts_pop0_joint = genotypes_pop0_joint.count_alleles()
	genotypes_pop1_joint = allel.HaplotypeArray(haps_pop1_joint).to_genotypes(ploidy=2)
	allele_counts_pop1_joint = genotypes_pop1_joint.count_alleles()

	sfs_joint = allel.joint_sfs(allele_counts_pop0_joint[:,1], allele_counts_pop1_joint[:,1])
	num_sites = sum(sum(sfs_joint))
	#print(ts.num_sites)
	sfs_joint = dadi.Spectrum(sfs_joint)
	sfs_joint.to_file(out_path)
	sfs_joint[0,0] = total_length - num_sites # need to get the number of nonvariant sites for the [0,0] entry
	sfs_joint.to_file(out_path_nonvariant)




def OoA_func(params, ns, pts):
	'''
	Out of Africa model for dadi to simulate SFS
	'''
	nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs = params

	xx = dadi.Numerics.default_grid(pts)

	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuB, m12=mAfB, m21=mAfB)

	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)

	nuEu_func = lambda t: nuEu0 * (nuEu/nuEu0)**(t/TEuAs)
	nuAs_func = lambda t: nuAs0 * (nuAs/nuAs0)**(t/TEuAs)
	phi = dadi.Integration.three_pops(phi, xx, TEuAs, nu1=nuAf, nu2=nuEu_func, nu3=nuAs_func,
                                    m12=mAfEu, m13=mAfAs, m21=mAfEu,
                                    m23=mEuAs, m31=mAfAs, m32=mEuAs)


	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx,xx))
	return fs

def compare_msprime_dadi_OutOfAfrica(input_fids, output_path, sample_size=20):
	#For parameter reference
	#p0 = [nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs]
	OoA_popt = [1.68, 0.287, 0.129, 3.74, 0.070, 7.29, 3.65, 0.44, 0.28, 1.40, 0.607, 0.396, 0.058]
	OoA_pts_l = [30, 40, 50]
	OoA_ns = [20,20,20]
	OoA_extrap_func = dadi.Numerics.make_extrap_func(OoA_func)
	OoA_model = OoA_extrap_func(OoA_popt, OoA_ns, OoA_pts_l)
	OoA_model = OoA_model.marginalize([2])

	msprime_joint_sfs = dadi.Spectrum([[0]*(sample_size+1)]*(sample_size+1))


	for fid in input_fids:
		msprime_joint_sfs_temp = dadi.Spectrum.from_file(fid)
		msprime_joint_sfs += msprime_joint_sfs_temp

	fig = plt.figure(219033)
	fig.clear()
	dadi.Plotting.plot_2d_comp_multinom(OoA_model, msprime_joint_sfs, vmin=1, resid_range=50, show=False)
	fig.savefig(output_path)



def OoA_2D_Af_Eu_func(params, ns, pts):
	'''
	Out of Africa model for dadi to simulate SFS
	'''
	nuAf, nuB, nuEu0, nuEu, mAfB, mAfEu, TAf, TB, TEuF = params

	xx = dadi.Numerics.default_grid(pts)

	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuB, m12=mAfB, m21=mAfB)

	nuEu_func = lambda t: nuEu0 * (nuEu/nuEu0)**(t/TEuF)
	phi = dadi.Integration.two_pops(phi, xx, TEuF, nu1=nuAf, nu2=nuEu_func, m12=mAfEu, m21=mAfEu)

	fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
	return fs


def fit_dadi_model(sfs_files,output_pdf_name,output_text_name,demo_model,fit_seed, sample_size=20):
	np.random.seed(int(fit_seed))

	msprime_joint_sfs = dadi.Spectrum([[0]*(sample_size+1)]*(sample_size+1))


	for fid in sfs_files:
		msprime_joint_sfs += dadi.Spectrum.from_file(fid)

	pts_l = [40, 30, 50]
	ns = msprime_joint_sfs.sample_sizes

	if demo_model == 'IM':
		func = dadi.Demographics2D.IM
		params = [0.5 , 2.0 , 2.0 , .5, 1, 1]
		lower_bound = [0, 0, 0, 0, 0, 0]
		upper_bound = [1, 50, 50, 10, 10, 10]

	if demo_model == 'IM_pre':
		func = dadi.Demographics2D.IM_pre
		params = [2, 0.5, 0.5 , 3, 3, .5, 1, 1]
		lower_bound = [0, 0, 0, 0, 0, 0, 0, 0]
		upper_bound = [50, 10, 1, 50, 50, 10, 10, 10]

	if demo_model == 'split_mig':
		func = dadi.Demographics2D.split_mig
		params = [1.5, 1.5, 0.5, 0.5]
		lower_bound = [0, 0, 0, 0]
		upper_bound = [50, 50, 10, 10]

	if demo_model == 'IM_fsc':
		func = dadi.Demographics2D.IM_fsc
		params = [1.5, 1.5, 0.5, 0.5, 0.5]
		lower_bound = [0, 0, 0, 0, 0]
		upper_bound = [50, 50, 10, 10, 10]

	if demo_model == 'Gute2pop':
		func = OoA_2D_Af_Eu_func
		params = [1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5]
		lower_bound = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		upper_bound = [50, 50, 50, 50, 10, 10, 10, 10, 10]


	extrap_function = dadi.Numerics.make_extrap_log_func(func)

	#Not useful right now
	try:
		fixed
	except:
		fixed = [None]*len(params)


	print('staring_optimization',fit_seed)

	p_guess = dadi.Misc.perturb_params(params, lower_bound=lower_bound, upper_bound=upper_bound)

	popt = dadi.Inference.optimize(p_guess, msprime_joint_sfs, extrap_function, pts_l, multinom=True, verbose = True,
					lower_bound=lower_bound, upper_bound=upper_bound,
					fixed_params=fixed)

	# #For debugging
	# popt = p_guess

	model = extrap_function(popt, ns, pts_l)

	ll = dadi.Inference.ll_multinom(model,msprime_joint_sfs)
	print('ll (seed '+str(fit_seed)+'):' + str(ll))

	theta0 = dadi.Inference.optimal_sfs_scaling(model,msprime_joint_sfs)
	print('theta0 (seed '+str(fit_seed)+'):'+str(theta0))


	fid = open(output_text_name,'w')
	fid.write('ll: ' + str(ll) + '\ntheta0: ' + str(theta0) + '\nParams:\n['+','.join([str(ele) for ele in popt])+']')
	fid.close()

	import matplotlib.pyplot as plt
	fig = plt.figure(219033)
	fig.clear()
	dadi.Plotting.plot_2d_comp_multinom(model, msprime_joint_sfs, vmin=1, show=False)
	fig.savefig(output_pdf_name)


def get_dadi_output_IM(indir,model,dadi_seeds,ofile):

	output = []
	for seed in dadi_seeds:
		path=indir+model+"/model_params_"+str(seed)+".txt"
		with open(path) as infile:
			file = infile.readlines()
			#print(file)
			ll = file[0].split()[1]
			theta = file[1].split()[1]
			nu1 = file[3].split(",")[0]
			nu1=nu1.replace("[","")
			nu2 = file[3].split(",")[1]
			T = file[3].split(",")[2]
			m1 = file[3].split(",")[3]
			m2 = file[3].split(",")[4]
			m2=m2.replace("]","")
			line = [ll,theta,nu1,nu2,T,m1,m2]
			output.append(line)

	output.sort(key=lambda x: float(x[0]), reverse=True)

	with open(indir+"/"+ofile, 'w') as outfile:
		outfile.write("ll\ttheta\tnu1\tnu2\tT\tm1\tm2\n")
		outfile.writelines('\t'.join(i) + '\n' for i in output)


def get_dadi_output_Gute2pop(indir,model,dadi_seeds,ofile):

	output = []
	for seed in dadi_seeds:
		path=indir+model+"/model_params_"+str(seed)+".txt"
		with open(path) as infile:
			file = infile.readlines()

			ll = file[0].split()[1]
			theta = file[1].split()[1]
			nuAf = file[3].split(",")[0]
			nuAf=nuAf.replace("[","")
			nuB = file[3].split(",")[1]
			nuEu0 = file[3].split(",")[2]
			nuEu = file[3].split(",")[3]
			mAfB = file[3].split(",")[4]
			mAfEu = file[3].split(",")[5]
			TAf = file[3].split(",")[6]
			TB = file[3].split(",")[7]
			TEuF = file[3].split(",")[8]
			TEuF=TEuF.replace("]","")
			line = [ll,theta,nuAf,nuB,nuEu0,nuEu,mAfB,mAfEu,TAf,TB,TEuF]
			output.append(line)

	output.sort(key=lambda x: float(x[0]), reverse=True)

	with open(indir+"/"+ofile, 'w') as outfile:
		outfile.write("ll\ttheta\tnuAf\tnuB\tnuEu0\tnuEu\tmAfB\tmAfEu\tTAf\tTB\tTEuF\n")
		outfile.writelines('\t'.join(i) + '\n' for i in output)


def get_best_dadi_runs(indir,seeds,outfile):

	header_count = 0
	with open(outfile, 'w') as ofile:
		for seed in seeds:
			with open(indir+"/"+str(seed)+"/dadi_analysis/dadi_results_sorted.txt") as infile:
				file = infile.readlines()
				if header_count == 0:
					ofile.write(file[0])
					header_count =+1
				ofile.write(file[1])
