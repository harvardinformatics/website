Title: RNA-Seq differential expression workshop
Date: 2016-03-22
Author: Tim Sackton
Category: Tutorials
Tags: RNA-Seq, Workshop, Sleuth, kallisto, DESeq2
Summary: This tutorial provides a workflow for RNA-Seq differential expression analysis using DESeq2, kallisto, and Sleuth

Setup
------

To run this workshop you will need:

1. R (https://cran.r-project.org/)
1. the DESeq2 bioconductor package (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
3. kallisto (https://pachterlab.github.io/kallisto/)
4. sleuth (pachterlab.github.io/sleuth/)

Introduction
------

There are two main approaches for detecting differential expression of genes and transcripts using RNA-seq data. The first, older, approach is based on first mapping reads to transcripts (using tools such as RSEM or Cufflinks), and then using the estimated counts of reads that map to each transcript or gene as the input to a statistical model, typically a negative binomial model of read counts, such as is implemented in R packages like edgeR or DESeq2. This approach has the advantage of very good accuracy, but is quite slow.

An alternate approach is based on the idea of pseudo-alignment or pseudo-mapping. In this case, reads are not directly aligned but rather k-mers are used to assign reads to the most likely transcript that generated them. This approach is incredibly fast as it does not have to do the time consuming computation of alignment statistics, and is nearly as accurate as gold-standard mapping approachs such as RSEM. kallisto and salmon are two software programs that use this approach, and sleuth is the main analysis package to analyze the data after generating transcript counts.

Today, we are going to work through a kallisto & sleuth example, and then compare it to the same data analyzed with RSEM and DESeq2. Because RSEM is slow, we will use precomputed count data, but we will run kallisto during the workshop.


Data
----

We are going to use as an example some RNA-seq data generated from male and female ostriches from two different tissues for this paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3603317/

This is a nice example data set, as it has a simple, balanced, two-factor experiment design, with three biological replicates per sex*tissue combination.

The fastq files and ostrich annotations necessary for this analysis are already on the cluster.


kallisto
---------

We will run kallisto in an interactive node on the Odyssey cluster, although for normal use you would probably want to run it as a batch job. If you do not have an Odyssey cluster account, please partner up with someone who does.

First we'll need to get an interactive node for this work:
	:::bash
	$ srun -p interact --pty --mem 4000 -t 0-1:00 /bin/bash

Then we need to load the modules on the cluster:

	:::bash
	$ source new-modules.sh
	$ module load gcc/4.8.2-fasrc01 openmpi/1.10.0-fasrc01 kallisto/0.42.4-fasrc01

Now we can run kallisto in two steps. The first is to build an index for the species we are going to analyze, in this case ostrich. To do this all you need is a fasta file of transcripts, such as what can be downloaded from UCSC, Ensembl, or NCBI.

We'll use the kallisto index command like so:
	:::bash
	$ kallisto index -i ostrich /n/regal/informatics/workshops/data/ostrich/ostrich.transcripts.fa 
	
This should only take a few minutes. In the mean time, does anyone have any questions so far?

Now we need to quantify our transcripts. For this example we will only run this for one sample (a single replicate of one tissue by sex combination), but for real analysis you'd be launching several different kallisto commands, one for each sample you want to quantify.

	:::bash
	$ kallisto quant -i ostrich -o BF1 -t 1 -b 1 /n/regal/informatics/workshops/data/ostrich/BF1.fq.1 /n/regal/informatics/workshops/data/ostrich/BF1.fq.2

This will take around 15 minutes, although doing this for a real analysis would likely take a bit longer. Partly this is because I am not doing *bootstrapping* (randomly resampling reads to estimate sampling variance), which is controlled by the -b flag. In real data, you'd likely want to do at least 100 bootstraps. 

While this runs, I will give a brief overview of the RSEM pipeline (read alignment) and discuss some of the issues associated with read counting. We will then turn to analyzing pre-generated RSEM and kallisto output for all 12 samples in R, using DESeq2 and sleuth.

Differential expression with DESeq2
------

First, we need to load the libraries we'll use.

	:::R
	library('DESeq2')
	library('RColorBrewer')

Now, we need to read our count data into a data frame so that we can analyze it with DESeq2. We can do this directly from web:

	:::R
	filePath <- "http://software.rc.fas.harvard.edu/ngsdata/workshops/RNAseq_2016/RSEM/"
	sampleNames <- c("BF1", "BF2", "BF3", "BM1", "BM2", "BM3", "LF1", "LF2", "LF3", "LM1", "LM2", "LM3")
	countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, ".isoforms.results"), header=T, sep="\t"), simplify=F)

countData.list now is a list of data frames, where each data frame contains the RSEM results from one sample. We need to select just the expected_count columns and rename them to be the sample names. We'll do this in a couple of steps.

	:::R
	countData.df <- do.call("cbind", countData.list)
	colsToKeep <- c(1,grep("expected_count", names(countData.df)))
	ct <- countData.df[,colsToKeep]
	names(ct) <- c("transcript_id", sampleNames)
	ct[,2:13] <- round(ct[,2:13])

Now we will construct a DESeq data object using the function DESeqDataSetFromMatrix:

	:::R
	sampleMetaData <- data.frame(tissue=c(rep(c("brain"), 6), rep(c("liver"), 6)), sex=rep(c(rep(c("F"), 3), rep(c("M"), 3)),2))
	rownames(sampleMetaData) <- sampleNames
	rsem.in <- DESeqDataSetFromMatrix(countData = ct, colData = sampleMetaData, design = ~ tissue + sex, tidy = T)
	rsem.de <- DESeq(rsem.in)

Look at how the dispersion is estimated:
	:::R
	plotDispEsts(rsem.de)

Start with some quality control: let's make a PCA plot of the samples

	:::R
	rsem.rlog <- rlog(rsem.de)
	plotPCA(rsem.rlog)
	plotPCA(rsem.rlog, intgroup=c("sex", "tissue"))
	
Given this, we might want to either analyze each tissue separately, or include an interaction term in our model. For simplicity we will focus on the brain samples.

	:::R
	rsem.brain <- DESeqDataSetFromMatrix(countData = ct[,1:7], colData = sampleMetaData[1:6,], design = ~ sex, tidy = T)
	rsem.brain.de <- DESeq(rsem.brain)
	rsem.brain.res <- results(rsem.brain.de)
	summary(rsem.brain.res)
	plotMA(rsem.brain.res, ylim=c(-3.5,1.5))
	plotCounts(rsem.brain, gene=which.min(rsem.brain.res$padj), intgroup="sex", pch=19)
	
There is a lot more that can be done with DESeq2, including varying a number of options from the defaults, but this should give you a brief introduction to the process of going from count tables to differentially expressed genes.

Differential expression with sleuth
------

First, load sleuth:

	:::R
	library(sleuth)
	
Now we need to download the kallisto data, which we will do using the command line:

	:::bash
	wget https://software.rc.fas.harvard.edu/ngsdata/workshops/kallisto_out.tar
	tar xf kallisto_out.tar 
	
Back in R:

	:::R
	base_dir <- "~/Projects/RC/workshops/2016-03-22-RNAseq/kallisto_out"
	sample_id <- dir(file.path(base_dir))
	kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
	sleuth.sampledata <- data.frame(sample=sample_id, tissue=sampleMetaData$tissue, sex=sampleMetaData$sex, path=kal_dirs, stringsAsFactors=F)
	
Run sleuth:

	:::R
	sleuth.all <- sleuth_prep(sleuth.sampledata, ~ tissue + sex + tissue:sex)
	sleuth.all <- sleuth_fit(sleuth.all)
	sleuth.all <- sleuth_wt(sleuth.all, 'sexM')
	sleuth.all <- sleuth_wt(sleuth.all, 'tissueliver')
	sleuth.all <- sleuth_wt(sleuth.all, 'tissueliver:sexM')	

Plot results:
	:::R
	plot_pca(sleuth.all, color_by="tissue", text_labels=T)
	
Similar to the RSEM results, there is a huge effect of tissue. Let's again focus on the brain data.

	:::R
	sleuth.brain.in <- sleuth.sampledata[1:6,]
	sleuth.brain <- sleuth_prep(sleuth.brain.in, ~ sex)
	sleuth.brain <- sleuth_fit(sleuth.brain)
	sleuth.brain <- sleuth_wt(sleuth.brain, 'sexM')
	sleuth.brain.res <- sleuth_results(sleuth.brain, 'sexM')

