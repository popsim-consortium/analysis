"""
Utilities for working with stairwayplot.
"""
import pathlib
import logging
import subprocess
import tempfile
import shutil
import concurrent.futures

import tskit
import allel
import numpy as np
import tqdm


def write_stairway_sfs(ts, sfs, path):
    """
    Writes the SFS to stairway plot format to the specified file.
    """
    # write stairwayplot input
    with open(path, "w") as out:
        # order is name,n_samples,sequence_length,lowest_sfs_bin,highest_sfs_bin
        print(
            "msp", ts.num_samples, int(ts.sequence_length), 1,
            ts.num_samples - 1, sep="\t", file=out)
        # TODO can we use numpy to do this more efficiently?
        for x in sfs:
            print(int(x), end="\t", file=out)
        print(file=out)


class StairwayPlotRunner(object):
    """
    Run stairway plots.
    """
    def __init__(self, workdir, stairway_dir, java_exe="java"):
        self.workdir = pathlib.Path(workdir)
        shutil.rmtree(self.workdir, ignore_errors=True)
        self.workdir.mkdir(parents=True)
        stairway_path = pathlib.Path(stairway_dir)
        self.classpath = "{}:{}".format(stairway_path, stairway_path / "swarmops.jar")
        self.java_exe = java_exe

    def ts_to_stairway(self, ts_path, num_bootstraps=1):
        """
        Converts the specified tskit tree sequence to text files used by
        stairway plot.
        """
        ts = tskit.load(ts_path)

        # count alleles, bootstrap over sites, return the SFS minus the 0% bin
        haps = np.array(ts.genotype_matrix())
        genotypes = allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
        allele_counts = genotypes.count_alleles()
        sfs = allel.sfs(allele_counts[:, 1])
        sfs = sfs[1:len(sfs)]

        # write stairwayplot input
        filename = self.workdir / "sfs_0.txt"
        write_stairway_sfs(ts, sfs, filename)
        stairway_files = [filename]

        # Write bootstrapped inputs
        for j in range(1, num_bootstraps + 1):
            nsites = np.shape(allele_counts)[0]
            bootset = np.random.choice(np.arange(0, nsites, 1), nsites, replace=True)
            bootac = allele_counts[bootset, :]
            bootsfs = allel.sfs(bootac[:, 1])
            bootsfs = bootsfs[1:len(bootsfs)]
            filename = self.workdir / "sfs_{}.txt".format(j)
            write_stairway_sfs(ts, bootsfs, filename)
            stairway_files.append(filename)
        return stairway_files

    def _run_theta_estimation(self, input_file):
        """
        Runs stairway plot on the specified file, resulting in the output being written
        to the file input_file + ".addTheta".
        """
        num_runs = 1
        dim_factor = 5000
        cmd = (
            f"{self.java_exe} -cp {self.classpath} Stairway_plot_theta_estimation02 "
            f"{input_file} {num_runs} {dim_factor}")
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)

    def run_theta_estimation(self, max_workers=None, show_progress=False):
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self._run_theta_estimation, infile)
                for infile in self.workdir.glob("sfs_*.txt")]
            with tqdm.tqdm(total=len(futures), disable=not show_progress) as progress:
                for future in concurrent.futures.as_completed(futures):
                    progress.update()

    def run_summary(self, output_file, mutation_rate, generation_time):
        """
        Runs stairway plot summary files in the work dir, writing the
        output to output_file and with the given parameters.
        """
        # First we need to create a temporary directory for the files, stairway expects
        # only these files to be present in the directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            for infile in self.workdir.glob("*.addTheta"):
                shutil.copy(infile, tmpdir)
            cmd = (
                f"{self.java_exe} -cp {self.classpath} "
                f"Stairway_plot_output_summary_commandline {tmpdir} "
                f"{mutation_rate} {generation_time} {output_file}")
            logging.info("Running:" + cmd)
            subprocess.run(cmd, shell=True, check=True)
