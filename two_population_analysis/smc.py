"""
Utilities for working with scm++
"""
import logging
import subprocess
import tskit


def write_smcpp_file(path):
    """
    Writes a smcpp input file given a treesequence
    """
    ts = tskit.load(path)
    # write a vcf intermediate input
    with open(path+".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file, 2)
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
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop1:"
    for n in range(ts.num_samples // 2):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)


def write_smcpp_file_twopops(path):
    """
    Writes a smcpp input file given a treesequence
    """
    ts = tskit.load(path)
    # write a vcf intermediate input
    with open(path+".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file, 2)
    # index the vcf
    # cmd = f"head -n100000 {path}.vcf > {path}small.vcf"
    # subprocess.run(cmd, shell=True, check=True)
    cmd = f"bgzip {path}.vcf"
    logging.info("Running:" + cmd)
    print (cmd)
    subprocess.run(cmd, shell=True, check=True)
    vz_file = f"{path}.vcf.gz"
    print("File = " + vz_file + "\n")
    cmd = f"tabix {vz_file}"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    # run smc++ for vcf conversion
    smc_file = f"{path}pop1.smc.gz"
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop1:"
    for n in range(ts.num_samples // 4):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    smc_file = f"{path}pop2.smc.gz"
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop2:"
    for n in range(ts.num_samples // 4, ts.num_samples // 2):
       cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    smc_file = f"{path}pop12.smc.gz"
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop1:"
    for n in range(ts.num_samples // 4):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    cmd = cmd + " pop2:"
    for n in range(ts.num_samples // 4, ts.num_samples // 2):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    smc_file = f"{path}pop21.smc.gz"
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop2:"
    for n in range(ts.num_samples // 4, ts.num_samples // 2):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    cmd = cmd + " pop1:"
    for n in range(ts.num_samples // 4):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)




"""
def run_smcpp_estimate(input_file, mutation_rate, ncores):
    cmd = (
        f"smc++ estimate "
        f"{mutation_rate} {input_file} --base {input_file} --cores {ncores}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
"""

def run_smcpp_estimate(input_file, base, mutation_rate, ncores):
    """
    Runs smc++ estimate on the specified file, resulting in the output being written
    to the file input_file.final.jason".
    """
    #InputFile = input_file
    cmd = (
        f"smc++ estimate "
        f"{mutation_rate} {input_file}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    #InputFile = input_file + "pop2.smc.gz"
    #cmd = (
    #       f"smc++ estimate "
    #       f" {mutation_rate} {InputFile}")
    #logging.info("Running:" + cmd)
    #subprocess.run(cmd, shell=True, check=True)


def run_smcpp_plot(input_file1, input_file2, output_file, generation_time):
    """
    Runs smc++ plot on the specified file, resulting in the output being written
    to the file input_file.png".
    """
    SplitFile = input_file1.split("/")
    SplitFile = SplitFile[:-1]
    SplitFile = "/".join(SplitFile)
    WorkingDir = SplitFile
    SplitDir = SplitFile + "/split"
    SplitFile = SplitFile + "/chr22.trees*.smc.gz"
    print ("The split file is " + SplitFile)
    #logging.info("Running:" + cmd)
    cmd = (
        f"smc++ split -o {SplitDir} {input_file1} {input_file2} {SplitFile}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = (
           f"smc++ plot {SplitDir}/joint.pdf {SplitDir}/model.final.json -c")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = (f"cp {SplitDir}/model.final.json {output_file}")
    subprocess.run(cmd, shell=True, check=True)

def get_smcpp_runs (indir,seeds,outfile):
    
    header_count = 0
    with open(outfile, 'w') as ofile:
        ofile.write("HarmonicMeanPop1\tHarmonicMeanPop2\tDivTime\n")
        for seed in seeds:
            with open(indir+"/"+str(seed)+"/split/joint.csv") as infile:
                RowNumber = 0
                HarmonicMeanPop2 = 1
                PastTime = 0
                CurTime = 0
                for line in infile:
                    line.rstrip()
                    print (line)
                    SplitLine = line.split(',')
                    print (SplitLine[0])
                    if (SplitLine[0] == "pop2"):
                        if RowNumber == 0:
                            PastTime = float(SplitLine[1])
                            PastNe = float(SplitLine[2])
                        else:
                            CurTime = float(SplitLine[1])
                            CurNe = float(SplitLine[2])
                            HarmonicMeanPop2 = HarmonicMeanPop2 + (CurTime - PastTime) / PastNe
                            print (str(HarmonicMeanPop2))
                            PastTime = CurTime
                            PastNe = CurNe
                        RowNumber = RowNumber + 1
                DivTime = PastTime
                HarmonicMeanPop2 = CurTime / HarmonicMeanPop2
            with open(indir+"/"+str(seed)+"/split/joint.csv") as infile:
                RowNumber = 0
                HarmonicMeanPop1 = 1
                EndFlag = 0
                for line in infile:
                    SplitLine = line.split(',')
                    if (SplitLine[0] == "pop1"):
                        if RowNumber == 0:
                            PastTime = float(SplitLine[1])
                            PastNe = float(SplitLine[2])
                        else:
                            CurTime = float(SplitLine[1])
                            CurNe = float(SplitLine[2])
                            if CurTime > DivTime:
                                HarmonicMeanPop1 = HarmonicMeanPop1 + (DivTime - PastTime) / PastNe
                                PastTime = CurTime
                                PastNe = CurNe
                                EndFlag = 1
                            else:
                                HarmonicMeanPop1 = HarmonicMeanPop1 + (CurTime - PastTime) / PastNe
                                PastTime = CurTime
                                PastNe = CurNe
                        RowNumber = RowNumber + 1
                    if EndFlag == 1:
                        break
            HarmonicMeanPop1 = DivTime/ HarmonicMeanPop1
            ofile.write(str(HarmonicMeanPop1) + "\t" + str(HarmonicMeanPop2) + "\t" + str(DivTime) + "\n")


