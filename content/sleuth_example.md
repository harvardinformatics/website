Title: Sleuth example on Odyssey
Date: 2015-12-10 00:00
Category: Blog
Tags: Next-Gen Sequencing Transcriptome RNA-seq Differential Expression Sleuth
Summary: An example of running a Sleuth analysis on Odyssey cluster


[Sleuth](http://pachterlab.github.io/sleuth/) is a tool for the analysis and comparison of multiple related RNA-Seq experiments for which transcript abundances have been quantified with [kallisto](http://pachterlab.github.io/kallisto/). It implements statistical algorithms for differential expression analysis, and provides tools for exploratory data analysis.


#### 1  Start SLURM session and R

Sleuth is an R package and can be installed in a user account on Odyssey cluster. We recommend running it in an interactive SLURM session on Odyssey.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):

	:::bash
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash

In the interactive session we can start R, after loading the appropriate module.
To use modules of the Lmod system (which is recommended), the following command is required:

	:::bash
	$ source new-modules.sh

Sleuth requires R version 3.2.1 or later, so we can load the following version:

	:::bash
	$ module load R/3.2.2-fasrc01 R_packages/3.2.2-fasrc03

In order for a user to install R packages in the user account environment the following R environment variable needs to be set:

	:::bash
	$ export R_LIBS_USER=$HOME/R:$R_LIBS_USER

Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs.
For example, if you create your working directory under regal and name it "sleuth_example", the full path to the directory could be: ~/regal/sleuth_example



#### 2  Download example data

Before starting R and installing sleuth, we need to get the example data. The example data are human fibroblast RNA-Seq data (full data available from GEO at accession GSE37704), six samples total, three biological replicates in each of two conditions (three samples of HOXA1 knockdown, marked "HOXA1KD", versus three control samples, marked "scramble"). These RNA-Seq expression results have already been quantified with kallisto.
 
Download the example data and uncompress in the working directory:

	:::bash
	$ wget http://bio.math.berkeley.edu/sleuth/cuffdiff2/cuffdiff2_data_kallisto_results.zip
	$ unzip cuffdiff2_data_kallisto_results.zip

This will produce: a directory called "results" with the kallisto quantification data for each of the six samples (in six sub-folders); and a table (textfile) "hiseq_info.txt" with information about the samples that sleuth needs.

Note the location of the folder where the data have been downloaded as this path will needed in the sleuth analysis.



#### 3  Install sleuth software

Now, we are ready to start R, by simply typing in the following:

	:::bash
	$ R

This will produce the R prompt.

The sleuth package can now be installed with the following. First install rhdf5 by typing:

	:::r
	source("http://bioconductor.org/biocLite.R")
	biocLite("rhdf5")

Then install devtools:

	:::r
	install.packages("devtools")


Install sleuth by typing:

	:::r
	devtools::install_github("pachterlab/sleuth")

Finally, load sleuth with:

	:::r
	library("sleuth")

Optionally, to verify versions of R, sleuth and other loaded software, run:

	:::r
	sessionInfo()

This example works with sleuth version 0.28.0.


#### 4  Sleuth analysis preliminaries

Before running the sleuth analysis proper, we need to define some intermediate variables to bind sample data to experimental design information.

First specify where the example data kallisto results are located. Set an R variable with the base directory where data were downloaded in section 2 above:

	:::r
	base_dir <- "~/regal/sleuth_example" ## susbtitute with your own path to the sleuth example directory

Get the list of sample IDs:

	:::r
	sample_ids <- dir(file.path(base_dir,"results"))


Get paths to the kallisto results folders associated with the sample IDs:

	:::r
	kallisto_dirs <- sapply(sample_ids, function(id) file.path(base_dir, "results", id, "kallisto"))


Next, load the sample information we downloaded above, "hiseq_info.txt", that describes the experimental design and the relationship between the kallisto directories and the samples, and create table that ties sample name with condition and path to kallisto results folder:

	:::r
	s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
	s2c <- dplyr::select(s2c, sample = run_accession, condition)
	s2c <- dplyr::mutate(s2c, path = kallisto_dirs)

View the resulting table with:

	:::r
	print(s2c)

You should a table like the following (with your own working directory path):

	     sample condition                                        path
	1 SRR493366  scramble ~/sleuth_example/results/SRR493366/kallisto
	2 SRR493367  scramble ~/sleuth_example/results/SRR493367/kallisto
	3 SRR493368  scramble ~/sleuth_example/results/SRR493368/kallisto
	4 SRR493369   HOXA1KD ~/sleuth_example/results/SRR493369/kallisto
	5 SRR493370   HOXA1KD ~/sleuth_example/results/SRR493370/kallisto
	6 SRR493371   HOXA1KD ~/sleuth_example/results/SRR493371/kallisto



#### 5  Run sleuth analysis proper

The “sleuth object” for the analysis can be now constructed. It requires three commands.

First, load the kallisto processed data into the object:

	:::r
	so <- sleuth_prep(s2c, ~ condition)

Then estimate  parameters for the sleuth response error measurement model:

	:::r
	so <- sleuth_fit(so)

And finally, perform differential analyis (testing):

	:::r
	so <- sleuth_wt(so, 'conditionscramble')

The three steps should take up to 6-7 minutes altogether.


To generate a table of analysis results:

	:::r
	results_table <- sleuth_results(so, 'conditionscramble')

The table contains information about differential expression of transcripts between the control condition "scramble" versus the experimental gene knockdown condition "HOXA1KD". 

The table is already sorted by inclreasing q-value (qval column); if it were not (or if we wanted to sort the table by a different column), we could sort the table using:

	:::r
	results_ordered <- results_table[order(results_table$qval),]

To check how many transcripts are differentially expressed between the two conditions at significance level of q-value less than or equeall to 0.05, type:

	:::r
	table(results_ordered$qval <= 0.05)

This should produce the following result, indicating that 7,708 transcripts are differentially expressed (DE) at q-value <= 0.05.

	FALSE  TRUE 
	43117  7708

To write out these 7,708 DE transcripts into a tab-delimited text file, do:

	:::r
	write.table( subset(results_ordered, qval <= 0.05), file='sleuth.DE_transcripts.qval_0.05.txt', sep="\t",row.names=F, quote=F)


Optionally, to view the results in a browser and do further exploratory data analysis, generate the Shiny webpage of the results with:

	:::r
	sleuth_live(so) 

The webpage has several menus with analysis options: for instance, check the "MA plot", "volcano plot", and "test table" (this has the same information as the "results_table" we generated above) under the "analyses" menu; the "PCA" and "sample heatmap" plots under the "maps" menu; and the "scatter plots" under the "diagnostics" menu.



#### 6 References

Lior Pachter (2015). [A sleuth for RNA-Seq](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/). Blogpost, August 17, 2015.



