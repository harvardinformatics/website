Title: Trinotate workflow example on Odyssey
Date: 2015-11-25 00:00
Category: Tutorials
Tags: Next-Gen Sequencing, Transcriptome, Trinotate
Summary: An example of running the Trinotate suite of tools to annotate a given transcriptome (in this case the mouse transcriptome).
Status: draft

[Trinotate](https://trinotate.github.io/) is a suite for the functional annotation of transcriptomes, particularly de novo assembled transcriptomes. It uses a number of different well referenced methods for functional annotation, including homology search against sequence databases (BLAST+/SwissProt), protein domain identification (HMMER/PFAM), and comparison to currently curated annotation databases (like eggNOG, and Gene Ontology terms).


#### 1  Software required

Download and uncompress Trinotate code:

	:::bash
	$ wget https://github.com/Trinotate/Trinotate/archive/v2.0.2.tar.gz -O Trinotate-2.0.2.tar.gz
	$ tar xvf Trinotate-2.0.2.tar.gz

You can mark the location of the Trinotate directory with an environment variable, for easy reference:  

	:::bash
	$ export TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2


Other software programs needed are already installed on the Odyssey cluster and can be loaded as modules.
Remember, if using modules of the Lmod system, that the following command is required:
	
	:::bash
	$ source new-modules.sh
	
Programs like NCBI BLAST+ and HMMER can be loaded like:

	:::bash
	$ module load blast/2.2.29+-fasrc01
	$ module load hmmer/3.1b1-fasrc01

SQLite (required for Trinotate database integration) is also already installed on the cluster.


[TransDecoder](http://transdecoder.github.io/), a program that finds coding regions within transcripts and returns the most likely longest-ORF peptide candidates, is also needed to preprocess the transcriptome data. This is already installed in the Lmod module system and can be loaded with:

	:::bash
	module load TransDecoder/2.0.1-fasrc01



#### 2  Databases required

Trinotate uses customized versions SwissProt and Pfam databases. These can be downloaded from the Trinotate website into the $TRINOTATE_HOME folder created above:

	:::bash
	$ wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_sprot.trinotate_v2.0.pep.gz
	$ wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz


To rename, uncompress and index the SwissProt database for Blast use:

	:::bash
	$ mv uniprot_sprot.trinotate_v2.0.pep.gz uniprot_sprot.trinotate.pep.gz
	$ gunzip uniprot_sprot.trinotate.pep.gz
	$ makeblastdb -in uniprot_sprot.trinotate.pep -dbtype prot

The Blast indexing should take about a minute.



To uncompress and prepare the Pfam database for use with 'hmmscan' do:

	:::bash
	$ gunzip Pfam-A.hmm.gz
	$ hmmpress Pfam-A.hmm

This should take 2-3 minutes.


Finally, retrieve and uncompress the latest Trinotate pre-generated resource SQLite database, which contains Uniprot(SwissProt and Uniref90)-related annotation information: 

	:::bash
	$ wget "https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz
	$ gunzip Trinotate.sqlite.gz

Place the Trinotate.sqlite database in the working directory. Create a separate working directory for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.




#### 3  Transcriptome data

These are the transcripts that we need to annotate. Normally a user would have obtained such a set from a data source, or as a result of their own analysis, e.g. a de novo assembled transcriptome from an RNA-Seq analysis.

In this example, we have used a mouse transcriptome from Ensembl. It can be downloaded and uncompressed with:

	:::bash
	$ wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
	$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz

The transcriptome data can be placed in a separate folder (the /n/regal file system is recommended), which we can mark with an environment variable, for easy reference:  

	:::bash
	$ export TRANS_DATA=/path/to/trasncriptome/data


In this example, we have somewhat simplified the sequence headers of the transcriptome fasta file and renamed it to "mouse38_cdna.fa". It contains 94,080 transcripts.


A tab-delimited text file indicating the Gene/Transcript relationships (format: "gene_id(tab)transcript_id"), is also needed for Trinotate. This file can be created using information from the transcriptome data source. In our example we obtained the table using sequence header annotation in the transcriptome fasta file. 
If you are using Trinity assemblies, you can generate this file using a Trinity script. 
We have named this file "gene_transcript_map.txt". 

The mouse example data can be copied from the Informatics reference genome directory on Odyssey: /n/regal/informatics_public/ref/ensembl/release-82/mus_musculus . This directory (/n/regal/informatics_public/ref) contains annotation datasets from other reference organisms that could potentially also be used as inputs for this Trinotate example workflow. 

The mouse example data can also be downloaded from [here](https://github.com/gmarnellos/Trinotate_example_supplement).

 

The transcriptome nucleotide sequences are processed with TransDecoder to produce the most likely longest-ORF peptide candidates. As this is a longer job we ran it as a SLURM job, here is the script:

	:::bash
	#!/bin/bash
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 5000
	#SBATCH -p serial_requeue
	#SBATCH -o transd.out
	#SBATCH -e transd.err
	#SBATCH -t 800

	source new-modules.sh
	module load TransDecoder/2.0.1-fasrc01

	TRANS_DATA=/path/to/trasncriptome/data

	TransDecoder.LongOrfs -t $TRANS_DATA/mouse38_cdna.fa

	TransDecoder.Predict -t $TRANS_DATA/mouse38_cdna.fa


The job took about seven hours. The main output file used downstream is "mouse38_cdna.fa.transdecoder.pep", containing 89,368 peptide sequences.

As there are tens of thousands of sequences in the transcriptome, we can break the nucleotide and peptide fasta files into volumes, using the [split_fasta.pl](https://github.com/gmarnellos/Trinotate_example_supplement/blob/master/split_fasta.pl) script.

	:::bash
	$ split_fasta.pl $TRANS_DATA/mouse38_cdna.fa  $TRANS_DATA/mouse38_cdna.vol  10000
	$ split_fasta.pl $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep  $TRANS_DATA/mouse38_pep.vol  10000




#### 4  Run the Blast and Pfam searches


It would be a good idea to run these searches in a separate working directory.

Blastx of transcript nucleotide sequences against SwissProt. Run as a SLURM job array for the transcript fasta volumes. The script is:

	:::bash
	#!/bin/bash
	#SBATCH -n 8
	#SBATCH -N 1
	#SBATCH --mem 5000
	#SBATCH -p serial_requeue
	#SBATCH -o bx_%A_%a.out
	#SBATCH -e bx_%A_%a.err
	#SBATCH -t 2000
	#SBATCH --array=1-10

	source new-modules.sh
	module load blast/2.2.29+-fasrc01 

	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2

	blastx -query $TRANS_DATA/mouse38_cdna.vol.${SLURM_ARRAY_TASK_ID}.fasta -db $TRINOTATE_HOME/uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.vol.${SLURM_ARRAY_TASK_ID}.outfmt6  


Blastx jobs in this SLURM array took from 3 hours to almost 24 hours.



Concatenate result volumes into a single results file:

	:::bash
	$ cat blastx.vol.*.outfmt6 > blastx.mouse.sprot.outfmt6



Blastp of transcript TransDecoder peptide sequences against SwissProt. Run as a SLURM job array for the fasta volumes. The script is:
 
	:::bash
	#!/bin/bash
	#SBATCH -n 8
	#SBATCH -N 1
	#SBATCH --mem 3000
	#SBATCH -p serial_requeue
	#SBATCH -o bp_%A_%a.out
	#SBATCH -e bp_%A_%a.err
	#SBATCH -t 800
	#SBATCH --array=1-9

	source new-modules.sh
	module load blast/2.2.29+-fasrc01 

	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2

	blastp -query $TRABS_DATA/mouse38_pep.vol.${SLURM_ARRAY_TASK_ID}.fasta  -db $TRINOTATE_HOME/uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.vol.${SLURM_ARRAY_TASK_ID}.outfmt6


Blastp jobs in this SLURM array took from 2 to 8 hours.


Concatenate result volumes into a single results file:

	:::bash
	$ cat blastp.vol.*.outfmt6 > blastp.mouse.sprot.outfmt6



Run HMMER against Pfam to identify protein domains; run as a SLURM job:

 
	:::bash
	#!/bin/bash
	#SBATCH -n 16
	#SBATCH -N 1
	#SBATCH --mem 25000
	#SBATCH -p general
	#SBATCH -o hmm.out
	#SBATCH -e hmm.err
	#SBATCH -t 2400

	source new-modules.sh
	module load hmmer/3.1b1-fasrc01

	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2

	hmmscan --cpu 16 --domtblout TrinotatePFAM.out $TRINOTATE_HOME/Pfam-A.hmm $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep > pfam.log


SLURM job took about 30 hours.



#### 5  OPTIONAL: Blast searches against Uniref90 and additional processing

Trinotate uses a customized version of the Uniref90 database. This can be downloaded from the Trinotate website into the $TRINOTATE_HOME folder created above:

	:::bash
	$ wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_uniref90.trinotate_v2.0.pep.gz

To rename, uncompress and index the Uniref90 database for Blast use:

	:::bash
	$ mv uniprot_uniref90.trinotate_v2.0.pep.gz uniprot_uniref90.trinotate.pep.gz
	$ gunzip uniprot_uniref90.trinotate.pep.gz
	$ makeblastdb -in uniprot_uniref90.trinotate.pep -dbtype prot

The Blast indexing should take about 30 minutes.


Break the transcriptome into volumes of 3,000 sequences to runs jobs in parallel, using the split_fasta.pl script (usage: split_fasta.pl <input_file> <output_volume_prefix> <num_sequences_per_volume>):

	:::bash
	$ split_fasta.pl $TRANS_DATA/mouse38_cdna.fa  $TRANS_DATA/mouse38_cdna.uniref90.vol  3000
	$ split_fasta.pl $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep  $TRANS_DATA/mouse38_pep.uniref90.vol  3000



Blastx of transcript nucleotide sequences against Uniref90. Run as a SLURM job array in partition "general" for the transcript fasta volumes. The script is:

	:::bash
	#!/bin/bash
	#SBATCH -n 32
	#SBATCH -N 1
	#SBATCH --mem 25000
	#SBATCH -p general
	#SBATCH -o bx_%A_%a.out
	#SBATCH -e bx_%A_%a.err
	#SBATCH -t 6000
	#SBATCH --array=1-32

	source new-modules.sh
	module load blast/2.2.29+-fasrc01 

	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2

	blastx -query $TRANS_DATA/mouse38_cdna.uniref90.vol.${SLURM_ARRAY_TASK_ID}.fasta -db $TRINOTATE_HOME/uniprot_uniref90.trinotate.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 > blastx.uniref90.vol.${SLURM_ARRAY_TASK_ID}.outfmt6  



Blastx jobs in this SLURM array took from about 10 hours to almost 30 hours (not inlcuding any time the jobs may have spent in submission queue, job restarts, etc.).



Concatenate result volumes into a single results file:

	:::bash
	$ cat blastx.uniref90.vol.*.outfmt6 > blastx.mouse.uniref90.outfmt6



Blastp of transcript TransDecoder peptide sequences against Uniref90. Run as a SLURM job array in partition "serial_requeue" for the fasta volumes. The script is:
 
	:::bash
	#!/bin/bash
	#SBATCH -n 32
	#SBATCH -N 1
	#SBATCH --mem 25000
	#SBATCH -p serial_requeue
	#SBATCH -o bp_%A_%a.out
	#SBATCH -e bp_%A_%a.err
	#SBATCH -t 800
	#SBATCH --array=1-30

	source new-modules.sh
	module load blast/2.2.29+-fasrc01 

	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2

	blastp -query $TRABS_DATA/mouse38_pep.uniref90.vol.${SLURM_ARRAY_TASK_ID}.fasta  -db $TRINOTATE_HOME/uniprot_uniref90.trinotate.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 > blastp.uniref90.vol.${SLURM_ARRAY_TASK_ID}.outfmt6


Blastp jobs in this SLURM array took from about 8 to 16 hours (not inlcuding any time the jobs may have spent in submission queue, job restarts, etc.).


Concatenate result volumes into a single results file:

	:::bash
	$ cat blastp.uniref90.vol.*.outfmt6 > blastp.mouse.uniref90.outfmt6


Other OPTIONAL processing includes running the following programs; they are already installed on Odyssey as Lmod or legacy modules: 

- SignalP to predict signal peptides (Lmod module signalp/4.1c-fasrc01)
- tmHMM to predict transmembrane regions (legacy module hpc/TMHMM2.0c)
- RNAMMER 1.2 to identify rRNA transcripts (Lmod module rnammer/1.2-fasrc01)

Below is an example SLURM script to run these three programs sequentially. Alternatively, they can be run in parallel in three separate SLURM scripts using the same SLURM parameters. Each program would need about 2-4 hours to finish the mouse example dataset (tmHMM takes the longest); altogether they would take 10-12 hours.

	:::bash
	#!/bin/bash
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 5000
	#SBATCH -p serial_requeue
	#SBATCH -o optional.out
	#SBATCH -e optional.err
	#SBATCH -t 600
	
	TRANS_DATA=/path/to/trasncriptome/data
	TRINOTATE_HOME=/path/to/trinotate/Trinotate-2.0.2
	
	source new-modules.sh
	

	## Run SignalP to predict signal peptides
	
	module load signalp/4.1c-fasrc01
	
	signalp -f short -n signalp.out $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep
	
	
	## Run RNAMMER 1.2 to identify rRNA transcripts
	
	module load rnammer/1.2-fasrc01
	
	$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl --transcriptome $TRANS_DATA/mouse38_cdna.fa --path_to_rnammer /n/sw/fasrcsw/apps/Core/rnammer/1.2-fasrc01/rnammer
	
	
	## Run tmHMM to predict transmembrane regions
	
	module load legacy
	module load hpc/TMHMM2.0c
	
	tmhmm --short < $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep  > tmhmm.out



#### 6  Combine search results and create annotation report

In the working directory combine the results of the Blast and Pfam searches above into the SQLite database.

Perl modules are needed for this processing. Load with:

	:::bash
	$ source new-modules.sh
	$ module load perl-modules/5.10.1-fasrc12


Begin populating the SQLite database by loading:

- Transcript nucleotide sequences 
- Peptide sequences (from TransDecoder)
- Gene/Transcript relationships ("gene_transcript_map.txt" file)

The command to initialize the database is:

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite init --gene_trans_map $TRANS_DATA/gene_transcript_map.txt --transcript_fasta $TRANS_DATA/mouse38_cdna.fa --transdecoder_pep $TRANS_DATA/mouse38_cdna.fa.transdecoder.pep


Load Blast results:

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_swissprot_blastp  blastp.mouse.sprot.outfmt6

	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_swissprot_blastx  blastx.mouse.sprot.outfmt6

OPTIONAL - Load Blast results from Uniref90 searches if you have them:

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_trembl_blastp  blastp.mouse.uniref90.outfmt6  
	
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_trembl_blastx  blastx.mouse.uniref90.outfmt6


Load Pfam results:

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_pfam  TrinotatePFAM.out


OPTIONAL - Load results of additional processing (tmHMM, SignalP, RNAMMER):

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_tmhmm  tmhmm.out
	
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_signalp  signalp.out
	
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_rnammer  mouse38_cdna.fa.rnammer.gff


Each of these database loading commands takes about a minute (or less) to run.


Create annotation report table; a parameter that can be specified is "-E" (maximum E-value for reporting best blast hit and associated annotations); the default value (0.00001) is stringent enough, so that is what is used here:

	:::bash
	$ $TRINOTATE_HOME/Trinotate Trinotate.sqlite report -E 0.00001 > trinotate_annotation_report.xls

When the optional results are not loaded, the report generation takes about 20-30 minutes; most of the mouse example dataset transcripts, close to 95%, have annotation in this report.

When the optional results are loaded as well, report generation takes about 3-4 hours; about 97% of the transcripts have annotation in this case.


The report table includes columns for:

- gene_id
- transcript_id
- sprot_Top_BLASTX_hit
- prot_id
- prot_coords
- sprot_Top_BLASTP_hit
- Pfam
- eggnog
- gene_ontology_blast
- gene_ontology_pfam


With optional results loaded, the following columns also have data in the report:

- TrEMBL_Top_BLASTX_hit
- RNAMMER (very few entries with mouse example dataset)
- SignalP
- TmHMM
- TrEMBL_Top_BLASTP_hit



Voilà!


#### 7 References

[Trinotate](https://trinotate.github.io/) has not been published as a paper but the last section of their webpage has a list of references for the individual software programs that Trinotate uses for functional annotation.
