"""
Utilities for working with fastsimcoal
"""
import dadi
import tskit


def dadi_to_fsc_sfs(sfs_files,dadi_out_path,fsc_out_path, sample_size=20):


    ## get data sfs files and make them into one joint sfs and save
    msprime_joint_sfs = dadi.Spectrum([[0]*(sample_size+1)]*(sample_size+1))

    for fid in sfs_files:
        msprime_joint_sfs += dadi.Spectrum.from_file(fid)

    msprime_joint_sfs.to_file(dadi_out_path)

    ## convert dadi 2D sfs to FSC 2D sfs
    ## NB: FSC joint format file names look like this: <prefix>_jointMAFpop1_0.obs
    ## Where the first pop specified is listed in the rows and the second pop
    ## specified is listed in the columns.
    with open(fsc_out_path, 'w') as outfile:
        outfile.write("1 observation\n")
        ## Format column headers (i.e. d0_0 d0_1 d0_2 .. d0_n for deme 0 up to sample size of n)
        outfile.write("\t" + "\t".join(["d0_"+str(x) for x in range(sample_size+1)]) + "\n")

        ## Format row headers
        row_headers = ["d1_" + str(x) for x in range(sample_size+1)]

        with open(dadi_out_path) as infile:
            ## Get the second line of the dadi-style sfs which contains the data
            row_data = infile.readlines()[1].split()
            row_size = sample_size + 1
            ## Slice the row data into evenly sized chunks based on the number of columns
            rows = [row_data[i:i + row_size] for i in range(0, len(row_data), row_size)]
            ## Write out each row to the file
            for i, row_head in enumerate(row_headers):
                outfile.write(row_head + "\t" + " ".join(rows[i]) + "\n")

def get_fsc_output(path_to_fsc_analysis,num_runs,ofile,fsc_model):

    output = []
    for run_num in range(1,num_runs+1):

        with open(path_to_fsc_analysis+"run"+str(run_num)+"/"+fsc_model+"/"+fsc_model+".bestlhoods") as infile:
            file = infile.readlines()
            header = file[0]
            output.append(file[1].split())
    # note that the column this is sorted on may need to be changed depending on how many parameters there are
    output.sort(key=lambda x: float(x[6]), reverse=True)

    with open(path_to_fsc_analysis+ofile, 'w') as outfile:
        outfile.write(header)
        outfile.writelines('\t'.join(i) + '\n' for i in output)

def get_best_fsc_runs(indir,seeds,outfile):

    header_count = 0
    with open(outfile, 'w') as ofile:
        for seed in seeds:
            with open(indir+"/"+str(seed)+"/fsc_analysis/fsc_results_sorted.txt") as infile:
                file = infile.readlines()
                if header_count == 0:
                    ofile.write(file[0])
                header_count =+1
                ofile.write(file[1])

def get_bestlhood_fsc(path_to_fsc_analysis,num_runs,fsc_model):
    
    output = []
    for run_num in range(1,num_runs+1):
        nrows = 0
        LineToKeep = ""
        with open(path_to_fsc_analysis+"run"+str(run_num)+"/"+fsc_model+"/"+fsc_model+".brent_lhoods") as infile:
            nrows = 0
            for line in infile:
                #print (line)
                SplitLine = line.split()
                ElementsInList = len (SplitLine)
                if (ElementsInList >= 7):
                    #print (SplitLine[7])
                    if (nrows == 1):
                        MaxLL = SplitLine[7]
                        LineToKeep = line
                        LineToKeep = SplitLine[1] + "\t" + SplitLine[2] + "\t" + SplitLine[3] + "\t" + SplitLine[4] + "\t" + SplitLine[5] + "\t" + SplitLine[6] + "\t" + SplitLine[7] + "\t" + SplitLine[7]
                    if (nrows > 1):
                        if (SplitLine[7] < MaxLL):
                            MaxLL = SplitLine[7]
                            LineToKeep = SplitLine[1] + "\t" + SplitLine[2] + "\t" + SplitLine[3] + "\t" + SplitLine[4] + "\t" + SplitLine[5] + "\t" + SplitLine[6] + "\t" + SplitLine[7] + "\t" + SplitLine[7]
                    nrows = nrows + 1
        OutputFile = open(path_to_fsc_analysis+"run"+str(run_num)+"/"+fsc_model+"/"+fsc_model+".bestlhoods","w")
        OutputFile.write("ANCSIZE\tNPOP1\tNPOP2\tTDIV\tMIG12\tMIG21\tMaxEstLhood\tMaxObsLhood\n")
        OutputFile.write(LineToKeep)
        OutputFile.close()


