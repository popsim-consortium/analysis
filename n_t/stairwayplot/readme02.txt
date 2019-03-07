Stairway plot package version 0.2

Release:
Nov. 3, 2014

Introduction:
This is a java program package for inferring detailed ancestral population size changes using the stairway plot method. The method first uses site frequency spectrum (SFS) from DNA sequence data to estimate a serial of population mutation rates (i.e. \theta=4*\mu*Ne, where \mu is the mutation rate per generation and Ne is the effective population size), assuming a flexible multi-epoch demographic model. Then it estimates the effective population size changes along the time based on the \theta estimations. Please see below for the usage descriptions for each program. We suggest the following steps for a typical analysis:
	(1) Create 200 or more input files for Stairway_plot_theta_estimation02.class. Each input file represents a SFS of a sample. For example, you can use the ms package (http://home.uchicago.edu/rhudson1/source/mksamples.html) to simulate DNA sequence samples with 200 replications. Then use ms_output_to_stairway_plot_input.class (see below) to convert the outputs into 200 SFSs for Stairway_plot_theta_estimation02.class, one for each replication. For real data, you can write the SFS to an input file (see input file format below) and then use Stairway_plot_input_bootstrap.class to create (e.g. 199) additional input files from bootstrap samples (the first bootstrap sample is always the original sample) assuming the independence between SNPs.
	(2) Run Stairway_plot_theta_estimation02.class with each SFS input and collect the output files (i.e. addTheta files). Users are encouraged to parallel this step on a computer cluster.
	(3) Run Stairway_plot_output_summary.class with the addTheta files collected in step (2), and output a summary of the estimated population size changes.
	(4) Use Stairway_plot_summary_plot.class or any plotting program you prefered to draw lines of the median, 2.5 and 97.5 percentiles of the population size measures (as Ys) against the corresponding time measure (as X), from the summary table outputted from step (3). 

The programs are provided AS IS without charge. Source codes are distributed under the RECEX SHARED SOURCE LICENSE. You can find a copy of the license with the package. 

Main programs:

1. Stairway_plot_theta_estimation02.class

This java program takes a input file presenting the SFS of a sample and output a serial of estimations of \theta assuming a multi-epoch model. It needs SwarmOps java library version 1.0 (or later) from the Hvass Laboratories, which can be downloaded from http://www.hvass-labs.org/projects/swarmops/java/. It also needs SFS_lnL02.class to be located within the same folder. It output a file with a name <input_file>.addTheta, which appends the \theta estimations to the input file, where <input_file> is name of the input file.

Usage:
Change directory to the folder where Stairway_plot_theta_estimation02.class, SFS_lnL02.class and the SwarmOps jar file (e.g. swarmops.jar) locate.

In Windows 
java -cp .;swarmops.jar Stairway_plot_theta_estimation02 <input_file> <numRuns> <dimFactor>

In Linux/Unix:
java -cp .:swarmops.jar Stairway_plot_theta_estimation02 <input_file> <numRuns> <dimFactor>

<input_file> is the name (including path, if not located in the current folder) of the input file.
<numRuns> specifies the number of runs for the optimization algorithm. The default value is 1. A larger number may potentially increase the accuracy of the estimation with a cost of longer estimation time.
<dimFactor> is the dimensional factor which determines the iterations of search per run, which equals to dimFactor*d, where d is the number of \theta to be estimated. The default value is 5000. A larger number may potentially increase the accuracy of the estimation with a cost of longer estimation time.

Input file format:
Columns are separated by TABs.
First row: the first 5 columns are mandatory
	1st col: population id or simulation task name or user's brief note
	2nd col: number of sequences in the sample (nseq)
	3rd col: length of sequence (L)
	4th col: the smallest size of SNP used for estimation. If all SNPs will be used, then this value is 1.
	5th col: the largest size of SNP used for estimation. If all SNPs will be used, then this value is nseq-1.
Second row: nseq-1 columns, each is a count of the SNPs of a given size (i.e. \xi). All counts need to be provided, even though some of them will not be used for \theta estimation.
	1st col: count of the SNPs of size 1 (i.e. \xi_1)
	2nd col: count of the SNPs of size 2 (i.e. \xi_2)
	...
	nseq-1_th col: count of the SNPs of size nseq-1 (i.e. \xi_(nseq-1))

Output file format:
Columns are separated by TABs.
The first two rows are the same as the input file (see above). 
Beginning from the third row are the records of intermediate results:
	2nd col: the number of groups of \theta estimated (ngroup)
	4th col: -log Likelihood
	beginning from the 5th col: groups of \theta estimated
	beginning from the ngroup+5_th col: the corresponding value of \theta per site estimated for each group
Beginning from the row started with "final model:" are the records for the final results: The number after the "final model:" is the -log Likelihood. Following that is one row presenting the groups of \theta estimated, and one row presenting the corresponding value of \theta (over the whole length L) estimated for each group.
 
2. Stairway_plot_output_summary.class

This java program takes output files (i.e. the addTheta files) from Stairway_plot_theta_estimation and summarizes the \theta estimations into the estimations of measures of time (in the expected number of mutation(s) per site or in years) and population size (in \theta per site or in individuals). It uses a GUI to ask for the output files from Stairway_plot_theta_estimation, assumed mutation rate per site per generation, assumed generation time (in years), and a output file name. 

Usage:
Change directory to the folder where Stairway_plot_output_summary.class and simpleGUI.jar locate.

In Windows 
java -cp .;simpleGUI.jar Stairway_plot_output_summary

In Linux/Unix:
java -cp .:simpleGUI.jar Stairway_plot_output_summary

Input file format: 
The same as the output file format of Stairway_plot_theta_estimation.class (i.e. the addTheta files).

Output file format:
Columns are separated by TABs.
1st row: the first row of the addTheta file (supposed to be the same for all addTheta files inputted)
2nd row: the number of addTheta files inputted
4th row: file names of the addTheta files
5th row: the corresponding final -log Likelihood from the addTheta files
6th row: title row of the summary table
Beginning from the 7th row are the estimated measures of time and population sizes. Every two rows represent a "step" of the stairway plot.
	1st col: time measured in the expected number of mutation(s) per site
	2nd col: corresponding \theta
	3rd col: the median of the population size measured in \theta per site
	4th col: the 2.5 percentile estimation of the population size measured in \theta per site
	5th col: the 97.5 percentile estimation of the population size measured in \theta per site
	6th col: time measured in years
	7th col: the median of the population size measured in individuals
	8th col: the 2.5 percentile estimation of the population size measured in individuals
	9th col: the 97.5 percentile estimation of the population size measured in individuals
 

Accessory programs:

1. ms_output_to_stairway_plot_input.class

This java program converts the output of the ms program to multiple input files for Stairway_plot_theta_estimation02.class. It uses a GUI to ask for the output file from ms, an ID of the population simulated (simulation task name or user's brief note), the length of each sequence simulated, the number of replications of simulation, and a file name (stem) for output.

Usage:
Change directory to the folder where ms_output_to_stairway_plot_input.class and simpleGUI.jar locate.

In Windows 
java -cp .;simpleGUI.jar ms_output_to_stairway_plot_input

In Linux/Unix:
java -cp .:simpleGUI.jar ms_output_to_stairway_plot_input

2. Stairway_plot_input_bootstrap.class

This java program creates multiple input files for Stairway_plot_theta_estimation02.class, based on bootstrap samples from one input file. It uses a GUI to ask for the original input file for Stairway_plot_theta_estimation02.class, and a number for additional input files to create. The default number is 199, with which the program will create a total of 200 input files with extensions .bootstrap<1-200>. The bootstrap1 file is the same as the original input file, and bootstrap<2-200> are based on bootstrap sampling. This program use the Binomial class from the Java library called Colt, which can be downloaded from http://acs.lbl.gov/ACSSoftware/colt/. Please note the following copyright notice from Colt: 
_____________________________________________________________________
Copyright (c) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose is hereby granted without fee, provided that the above copyright notice appear in all copies and that both that copyright notice and this permission notice appear in supporting documentation. CERN makes no representations about the suitability of this software for any purpose. It is provided "as is" without expressed or implied warranty.
_____________________________________________________________________

Usage:
Change directory to the folder where Stairway_plot_input_bootstrap.class, simpleGUI.jar and colt.jar locate.

In Windows 
java -cp .;simpleGUI.jar;colt.jar Stairway_plot_input_bootstrap

In Linux/Unix:
java -cp .:simpleGUI.jar:colt.jar Stairway_plot_input_bootstrap

3. Stairway_plot_summary_plot.class

This java program plots the output of Stairway_plot_output_summary.class. It needs GRAL java library version 0.1 (or later) from the GRAL project, which can be downloaded from http://trac.erichseifert.de/gral/wiki/Download. It plots Ne (in 1k individuals) against Time (in 1k years). Users have the options to change some plotting settings, including the min and max for the x-axis (in 1k years) and y-axis (in 1k individuals), tick spacing for the x-axis and y-axis, and font size. After the plot is drawn on the screen, users can save the plot into files by accessing the mouse-right-click menu.

Usage:
Change directory to the folder where Stairway_plot_summary_plot.class, simpleGUI.jar and the GRAL jar file (e.g. gral-core-0.10.jar and VectorGraphics2D-0.9.1.jar) locate.

In Windows 
java -cp .;simpleGUI.jar;gral-core-0.10.jar;VectorGraphics2D-0.9.1.jar Stairway_plot_summary_plot

In Linux/Unix:
java -cp .:simpleGUI.jar:gral-core-0.10.jar:VectorGraphics2D-0.9.1.jar Stairway_plot_summary_plot

Test run:

This an example of test run on a PC with a Windows system
1. Unzip all files in stairway_plot_v0.2.zip to a folder, e.g. Stairway_plot. Unzip all files in testdata.zip to folder Stairway_plot\testdata.
2. Enter command line environment, change folder to Stairway_plot.
3. Run batch02.bat. 
4. Run "java -cp .;simpleGUI.jar Stairway_plot_output_summary" and choose all the .Theta files in Stairway_plot\testdata. Output results to Stairway_plot\two-epoch.summary.
5. Run "java -cp .;simpleGUI.jar;gral-core-0.10.jar;VectorGraphics2D-0.9.1.jar Stairway_plot_summary_plot" and choose Stairway_plot\two-epoch.summary.

Contact:
Xiaoming Liu, Ph.D.
Assistant Professor,
Human Genetics Center,
School of Public Health,
The University of Texas Health Science Center at Houston.
Email: xmliu.uth{at}gmail.com 

Citation:
Liu X and Fu YX. (2015) Exploring population size changes using SNP frequency spectra. Nature Genetics. 47(5):555-559.  
If you'd like to use the package to analyze your data, please contact the author (Xiaoming Liu).
