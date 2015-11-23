Title: Scripting Tips for Bioinformatics
Date: 2014-02-21 00:00
Category: Blog
Tags: Linux
Summary: Linux command line scripting tips  

With the advent of NextGen sequencing, many fields of bioinformatics have moved into the realm of high-performance computing. Some algorithms require hundreds of gigabytes of RAM, while others may take weeks to complete. Meanwhile the raw intensity files of a few sequencing runs can consume terabytes of storage. In this situation, its good to have a few scripting tricks up your sleeve. 


### Use ls -lh to view files sizes in human-friendly format
It can be difficult to read file sizes when they are 10+ digits. The -h option will display the number of bytes with an easy-to-read suffix such as "G" for gigabytes or "M" for megabytes. 

### Track memory usage with du -h and df -h
The du and df functions show you the memory usage of the current directory, and the current disk, respectively. The -h option displays memory counts in a human-friendly format. Use these functions to avoid being taken by surprise by "over quota" or "disk full" errors. 

### Use grep to search a file without loading it into RAM
When log files are too large to scan by eye, use grep to search for errors and warnings:

    :::shell-session
    $ grep -C2 -Ein 'error|exception|warning' log.txt

In addition, this command will print a few lines of text around the line where the error was found. This will help you determine the context in which the error occurred. 

### Use head and tail to view the beginning and end of files
For example, these commands display the first 50 lines and last 25 lines of mydata.txt:

    :::shell-session
    $ head -n 50 mydata.txt
    $ tail -n 25 mydata.txt

Like grep, head and tail allow you to view segments of files without loading the whole file into memory. 

### Use the find command to prune output files
Many bioinformatics algorithms generate intermediate data. Some generate thousands of small files, while others generate a few large ones. You can save space by deleting these files and retaining only the output you can use. For example, if you ran Trinity to assemble a transcriptome and wanted to delete all generated files except for the transcripts in Trinity.fasta and the runtimes in Trinity.timing, you could use the following command:

    :::shell-session
    $ find . -type f -not -name "Trinity.*" -print | xargs -I{} rm -vfr {}

### Use gunzip -c to examine compressed data without decompressing the entire file
Suppose you wanted to check whether a compressed FASTQ file used Phred base 33 or 64 quality scores. You could examine the first few lines of quality scores via:

    :::shell-session
    $ gunzip -c myreads.fastq.gz | head

After the first few lines are output, head sends a kill signal to the gunzip process, which stops decompressing the file. 

### Use the tar function "in flight" to move many files efficiently
The data transfer rate when moving thousands of small files can be slow as permission checks are performed for each file. The transfer rate can be dramatically increased by tarring the directory "in flight" as you copy it to its new location. An added advantage is that there's no need to store a tarred copy of the data in the original location:

    :::shell-session
    $ tar -C local_source_dir -cpf - . | tar -C /path/to/destination_dir -xvf -

### Use screen to avoid being disconnected from your server while running a command
An innocuous-looking call to find or du can take time to execute as the process generates sub-processes and traverses the directory structure. If you are concerned you might get disconnected, and thereby logged out, before your command has completed, run the command within a screen session. The session will persist even if you get disconnected from your server. Use

    :::shell-session
    $ screen -S <session_name>

to start a screen session, Ctrl-a Ctrl-d to disconnect from a session, and

    :::shell-session
    $ screen -r <session_name>

to rejoin a session you already created. To kill a session you are in, type Ctrl-D. To view a list of existing sessions, type:

    :::shell-session
    $ screen -ls