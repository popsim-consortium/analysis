"""
Utilities for working with scm++
"""
import logging
import subprocess
import tskit


def write_smcpp_file(path, mask_file=None):
    """
    Writes a smcpp input file given a treesequence
    """
    ts = tskit.load(path)
    chr_name = path.split(".")[0].split("/")[-1]
    # write a vcf intermediate input
    with open(path+".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file, 2, contig_id=chr_name)
    # index the vcf
    cmd = f"bgzip {path}.vcf"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    vz_file = f"{path}.vcf.gz"
    cmd = f"tabix {vz_file}"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    # run smc++ for vcf conversion
    smc_file = f"{path}.smc.gz"
    if mask_file:
        inter_mask_path = f"{path}.mask.gz"
        cmd = f"cat {mask_file} | bgzip > {inter_mask_path}"
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"tabix -p bed {inter_mask_path}"
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"smc++ vcf2smc --mask {inter_mask_path} {vz_file} "
        cmd = cmd + f"{smc_file} {chr_name} pop1:"
    else:
        cmd = f"smc++ vcf2smc {vz_file} {smc_file} {chr_name} pop1:"

    for n in range(ts.num_samples // 2):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)


def run_smcpp_estimate(input_file, base, mutation_rate, ncores):
    """
    Runs smc++ estimate on the specified file, resulting in the output being written
    to the file input_file.final.jason".
    """
    cmd = (
        f"smc++ estimate "
        f"--base {base} --cores {ncores} {mutation_rate} {input_file}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)


def run_smcpp_plot(input_file, output_file, generation_time):
    """
    Runs smc++ plot on the specified file, resulting in the output being written
    to the file input_file.png".
    """
    cmd = (
        f"smc++ plot {input_file}.png {input_file} -g {generation_time} -c")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = (f"cp {input_file}.csv  {output_file}")
    subprocess.run(cmd, shell=True, check=True)
