Title: BLAST and HMMER Databases on Odyssey
Date: 2019-04-22
Author: Nathan Weeks
Category: Software
Tags: Genome Databases, Odyssey
Summary: BLAST and HMMER reference databases on the Odyssey cluster

Reference databases for BLAST (nr, nt, uniref90) and HMMER (Pfam-A) are available on Odyssey.

## BLAST

The BLAST reference databases on Odyssey are formatted in version 5 to allow [limiting searches by taxonomy](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/blastdbv5.pdf), and require BLAST+ 2.8.0 or newer.

### Find / load BLAST software

Search for available NCBI BLAST software with the [module-query](https://www.rc.fas.harvard.edu/resources/documentation/software-on-odyssey/intro/) command:

    ```
    $ module-query ncbi-blast
    ...
        Versions:
            ncbi-blast/2.8.1+-fasrc01............... Core Built from the pre-compiled binary
            ncbi-blast/2.7.1+-fasrc01............... Core Built from the pre-compiled binary
            ncbi-blast/2.6.0+-fasrc01............... Core
            ncbi-blast/2.2.29+-fasrc03.............. Core Built for CentOS 7
    ...
    ```

Load the desired BLAST version:

    ```
    $ module load ncbi-blast/2.8.1+-fasrc01
    ```
     
### List / use available BLAST databases     
    
Loading the ncbi-blast module sets the BLASTDB environment variable to a directory containing available BLAST databases:
    
    ```
    $ blastdbcmd -list ${BLASTDB} -remove_redundant_dbs
    /n/holylfs/INTERNAL_REPOS/INFORMATICS/ref/blast/201907/nr Protein
    /n/holylfs/INTERNAL_REPOS/INFORMATICS/ref/blast/201907/nt Nucleotide
    /n/holylfs/INTERNAL_REPOS/INFORMATICS/ref/blast/201907/uniref90 Protein
    ```

To print just the database name to be used in conjunction with the BLAST -db option, and the database type:

    ```
    $ blastdbcmd -list ${BLASTDB} -remove_redundant_dbs | xargs -I {} basename {}
    nr Protein
    nt Nucleotide
    uniref90 Protein
    ```

To use one of the databases identified by the previous command, specify the `-db DB` option in the blast command in your job script:

    ```
	module load ncbi-blast/2.8.1+-fasrc01
    blastp -db nr ...
    ```

Multiple databases may be specified using a space-separated list as the `-db` option-argument (note this list must be surrounded by single or double quotes):

    ```
    blastp -db 'nr uniref90' ..
    ```

## HMMER

Search for / load availble [HMMER](http://www.hmmer.org/) versions:

    ```
    $ module-query hmmer
    ...
    hmmer : hmmer/3.1b2-fasrc01
    ...
    $ module load hmmer/3.1b2-fasrc01
    ```

List available HMMER profile databases:

    ```
    $ find db/2019Q2 -name '*.hmm.*' | sed -e 's#.*/##' -e 's#\.[^\.]*$##'  | uniq
    Pfam-A.hmm
    ```

Search a HMMER profile database with query sequences (TODO in the job script):

    ```
    module load hmmer/3.1b2-fasrc01
    hmmscan ${HMMERDB}/Pfam-A.hmm query.fa
    ```

## Database Update and Retention Policy

BLAST and HMMER databases are updated quarterly.
The two most-recent updates are retained, except directories with the "-stable" suffix will be retained for a longer time period.

Loading an ncbi-blast or hmmer environment module sets the BLASTDB environment variable and HMMERDB environment variables, respectively, to the most recent version:

    ```
    $ echo ${BLASTDB}
    /n/holylfs/INTERNAL_REPOS/INFORMATICS/ref/blast/201907
    ```

To see all available BLAST database versions, issue the following command (substitute the HMMERDB environment variable for HMMER):

    ```
    $ ls $(dirname ${BLASTDB})
    201801-stable
    201901-stable
    201904
    201907
    ```

To use a previous version, set the BLASTDB (HMMERDB) environment variable in your job script before invoking the appropriate BLAST (HMMER) command:

    ```
    export BLASTDB=$(dirname ${BLASTDB})/201801-stable
    blastp -db nr ...
    ```
