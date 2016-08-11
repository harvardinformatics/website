Title: De Novo Assembly: Trinity vs. Oases
Date: 2014-01-21 00:00
Category: Software
Tags: RNA-Seq Analysis,De Novo Assembly,Trinity,Oases
Summary: Oases and Trinity are two frequently used software packages for assembling short reads into transcripts. Both tools contain several complex steps and are difficult to evaluate on the basis of algorithm alone. For this reason, we compare the transcripts these tools create.

Oases and Trinity are two frequently used software packages for assembling short reads into transcripts. Both tools contain several complex steps and are difficult to evaluate on the basis of algorithm alone. For this reason, we compare the transcripts these tools create when run on the same set of Zebra finch, paired-end, 100 bp reads. 

First, Oases and Trinity were run on the cluster with slurm [scripts](https://github.com/harvardinformatics/bioblog/tree/master/Trinity.Oases.Assembly). Upon completion, each tool generates a fasta file (extension .fasta or .fa), which contains assembled transcripts. Using this perl and R [code](https://github.com/harvardinformatics/bioblog/tree/master/Summarize.Fasta), we can see some basic information about these assemblies in the plots below: 

![Oases.ZF_pUJ_n2_i5]({filename}/images/Oases.ZF_pUJ_n2_i5.png)  

The Oases assembly (above) contained three times as many transcripts, the Trinity assembly (below) contained fewer, longer transcripts. The difference in length distributions is reflected in the N50 values for the Trinity and Oases assemblies, which were 727 bp and 258 bp, respectively. 

![Trinity.ZF_pUJ_n2_i5]({filename}/images/Trinity.ZF_pUJ_n2_i5.png)

Transcripts from both tools were aligned to the reference genome with gmap, a genome aligner [3]. We then examined the extent to which known genes were covered by assembled transcripts, and the number of transcripts that did not map to any known genes. 

Roughly half of the 350,000 known exons in the Zebra finch genome were at least 30% covered by transcripts assembled by Oases and Trinity. However, as shown in the plot below, at higher coverage levels Trinity out-performed Oases in the number of known exons covered. Oases found approximately twice as many novel transcripts, which had no overlap with know exons. 

![Trinity.vs.Oases.Exon.Coverage.ZF_pUJ_n2_i5]({filename}/images/Trinity.vs_.Oases_.Exon_.Coverage.ZF_pUJ_n2_i5.png)


In conclusion, Trinity assembled fewer, longer transcripts which covered more of the reference genome. We may conclude that for analysis involving known genes, Trinity is the better tool for this data. However it may be that Oases is better for discovering novel transcripts, since the number of false positives among discovered transcripts is unknown.

1.  Schultz et al. [Oases: Robust de novo RNA-seq assembly across the dynamic range of expression levels](http://www.ncbi.nlm.nih.gov/pubmed/22368243). Bioinformatics. 2012 Feb 12.
2.  Haas et al. [De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis](http://www.nature.com/nprot/journal/v8/n8/full/nprot.2013.084.html). Nature Protocols. 2013 Aug;8(8):1494-512.
3.  Wu et al. [GMAP: a genomic mapping and alignment program for mRNA and EST sequences](http://bioinformatics.oupjournals.org/cgi/content/full/21/9/1859). Bioinformatics 2005, 21:1859-1875.