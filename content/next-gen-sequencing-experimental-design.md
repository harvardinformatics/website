Title: Next-Gen Sequencing Experimental Design
Date: 2015-01-01
Category: Tutorials
Author: Chris Williams
Tags: Next-Gen Sequencing
Summary: This tutorial provides useful tips for setting up sequencing experiments to provide the best results after analysis.

When planning a sequencing experiment, the number of reads, read type (single vs. paired-end) and number of replicates must be selected. Below are some considerations to bear in mind when making these decisions.

#### Sequencing Depth vs. Number of Replicates

The ability to detect transcripts in an RNA-Seq sample increases with sequencing depth. However, the advantage of additional depth decreases as more of a sample's fragments are covered. At some point, the vast majority of fragments will be covered at sufficient depth for detection, and the benefit of additional reads must be weighed against cost and the number of samples that can be sequenced with a given budget. 

A recent [study of bacterial genomes](http://www.biomedcentral.com/1471-2164/13/734) found that 5-10 million reads was sufficient to detect the majority of transcripts, while 2-3 million reads was sufficient to robustly detect differential expression. 

A [study of the chicken genome](http://www.biomedcentral.com/1471-2105/12/S10/S5) found that 30 million 75-bp reads were needed to detect all annotated chicken genes, while 10 million reads could detect about 80% of annotated genes. 

For expression quantification of mammalian genomes, 30-40 million reads may be sufficient, while lowly expressed isoforms may require [100-200 million 2 x 76 bp reads.](http://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf),

For the sake of thoroughness, it may be tempting to simply sequence 200 million reads for each sample, whatever the experimental question at hand. However, the experimental question may be better answered by spending the same amount of resources on sequencing more samples, each at less sequencing depth. 

[One study examined this tradeoff](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904521/) in an experiment that sought to identify differentially expressed genes between estradiol-treated and control breast cancer cells. They found that beyond 10 million reads per sample, the power to detect differentially expressed genes was greater when additional reads were used to sequencing more samples, rather than sequencing each sample at greater depth.

#### Distributing Samples Across Lanes

Since a single lane on an typical machine produces 200-400 million reads, multiple lanes may be needed to sequence all samples to the desired depth. For machines with a single lane, such as Illumina's NextSeq, multiple runs of the machine may be need. 

When using more than one lane of a single run for an experiment, it is a good idea to pool samples and run this same pool across multiple lanes, or multiple runs. There are several advantages to this approach. First, any batch effects that may stem from a single lane or run are avoided. This is especially true if one lane in the run fails and must be rerun at a later time. Second, if there are samples that are low-diversity, mixing the samples will increase the diversity and improve the machine's ability to focus the image of the flow-cell, identify clusters and call bases. If the cluster density of the pool is too low or too high, the concentration can be adjusted in future runs. An additional benefit of pooling of samples is the flexibility it lends to experiment design. After an initial data analysis, the sample pool can be sequenced at greater depth by running it again on additional lanes, while still avoiding batch effects.