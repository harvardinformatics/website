Title: Sequencing Core - MiniLIMS Billing Operations
Date: 2016-01-01
Category: Tutorials
Tags: Bauer Core, MiniLIMS
Author: Michele Clamp
Summary: How to for MiniLIMS sequencing billing

Sequencing billing is unavoidably complicated.   As well as the info provided by the submitters (Submissions and Samples) we need run info and samplesheet info pre-loaded by running scripts.    Christian will mark each Submission as BILLABLE as appropriate and then we can go to the InvoiceGenerator.

## Pre-Billing script 1 - the RunInfo files

Each run deposits a file with info about the name, the flowcell and other details.  Each of these needs to be uploaded.  As scanning the filesystem takes a while I generally just run it on directories in the current month.  The command looks like

    :::bash
    for i in /n/seqcfs/sequencing/analysis_finished/1501*; do 
        php plugins/Illumina/scripts/load_RunInfoFile.php -d $i -s ; 
    done

where 1501 is yymm and replace as necessary.    For testing purposes just remove the -s from the php line to check all is well.

## Pre-Billing script 2 - the SampleSheet files

These are very important as they contain which samples from which submissions were in each run.   These entries are then used to calculate how many lanes each submission occupied.

    :::bash
    for i in /n/seqcfs/sequencing/analysis_finished/1501*; do 
      php plugins/Illumina/NewSampleSheet.php -f $i/SampleSheet.csv -s
    done

As before 1501 is yymm and replace as necessary.    For testing purposes just remove the -s from the php line to check all is well.  This script can take a while to complete.