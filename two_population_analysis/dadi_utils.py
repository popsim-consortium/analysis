
import allel
import numpy as np
import msprime
from stdpopsim import homo_sapiens
import dadi
import os
import matplotlib.pyplot as plt
import random

def msprime_to_dadi_simulation_OutOfAfrica(path, seed, chrom, sample_size=20):
	'''
	Generate however many different SFS with msprime and convert+save them into SFS for dadi to use.
	'''
	#For testing
	# print(path, seed, chrom, sample_size)
	chrom = homo_sapiens.genome.chromosomes[chrom]
	model = homo_sapiens.GutenkunstThreePopOutOfAfrica()

	samples_pops_joint = [msprime.Sample(population=0, time=0)] * sample_size + [msprime.Sample(population=1, time=0)] * sample_size
	ts_pops_joint = msprime.simulate(
		samples=samples_pops_joint,
		recombination_map=chrom.recombination_map(),
		mutation_rate=chrom.default_mutation_rate,
		random_seed=seed,
		**model.asdict())
	haps_pops_joint = np.array(ts_pops_joint.genotype_matrix())

	#Break up the haplotypes into seperate populations based on sample_size
	haps_pop0_joint = haps_pops_joint[:,:sample_size]
	haps_pop1_joint = haps_pops_joint[:,sample_size:]

	genotypes_pop0_joint = allel.HaplotypeArray(haps_pop0_joint).to_genotypes(ploidy=2)
	allele_counts_pop0_joint = genotypes_pop0_joint.count_alleles()
	genotypes_pop1_joint = allel.HaplotypeArray(haps_pop1_joint).to_genotypes(ploidy=2)
	allele_counts_pop1_joint = genotypes_pop1_joint.count_alleles()

	sfs_joint = allel.joint_sfs(allele_counts_pop0_joint[:,1], allele_counts_pop1_joint[:,1])
	sfs_joint = dadi.Spectrum(sfs_joint)

	sfs_joint.to_file(path)


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

def compare_msprime_dadi_OutOfAfrica(input_fids, output_path):
	#For parameter reference
	#p0 = [nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs]
	OoA_popt = [1.68, 0.287, 0.129, 3.74, 0.070, 7.29, 3.65, 0.44, 0.28, 1.40, 0.607, 0.396, 0.058]
	OoA_pts_l = [30, 40, 50]
	OoA_ns = [20,20,20]
	OoA_extrap_func = dadi.Numerics.make_extrap_func(OoA_func)
	OoA_model = OoA_extrap_func(OoA_popt, OoA_ns, OoA_pts_l)
	OoA_model = OoA_model.marginalize([2])

	msprime_joint_sfs = dadi.Spectrum([[0]*21]*21) 

	for fid in input_fids:
		msprime_joint_sfs_temp = dadi.Spectrum.from_file(fid)
		msprime_joint_sfs += msprime_joint_sfs_temp

	fig = plt.figure(219033)
	fig.clear()
	dadi.Plotting.plot_2d_comp_multinom(OoA_model, msprime_joint_sfs, vmin=1, resid_range=50, show=False)
	fig.savefig(output_path)


def fit_dadi_model(sfs_files,output_pdf_name,output_text_name,demo_model,fit_seed):
	np.random.seed(int(fit_seed))

	msprime_joint_sfs = dadi.Spectrum([[0]*21]*21) 

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
					fixed_params=fixed,
					maxiter=2)

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