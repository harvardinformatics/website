Title: Intermediate Unix Workshop
Date: 2015-10-13 13:51
Category: Tutorials
Author: mclamp@harvard.edu
Tags: Linux, Workshop
Summary: This tutorial is intended for people who are familiar with the basics of unix but want to learn more about manipulating files and running commands.

[TOC]


## 1. Introduction and prerequisites

This tutorial is intended for people who are familiar with the basics of unix but want to learn more
about manipulating files and running commands.

*   A Harvard FAS RC cluster account 
    [https://account.rc.fas.harvard.edu/request/](account-request>)

*   A terminal program so you can log into the cluster 
    [https://rc.fas.harvard.edu/resources/access-and-login/](access-and-login>)

*   The ability to log into the cluster 
    [https://rc.fas.harvard.edu/resources/access-and-login/](access-and-login>)

*   <b>Familiarity with basic unix commands involving files (cp, rm , mv, cat, less, ls)
and/or have attended the ‘Basic Unix’ workshop previously.</b>




## 2. Summary of Commands Covered

* scp
* redirection (> and 2>&1)
* sort
* wc
* pipes |
* grep
* find
* awk
* bash for loops

### Functionality Covered

We call this the intermediate tutorial and it is intended for people who are comfortable with creating
files and moving around the filesystem but want to move onto more manipulation of data files.

## 3.  Transferring Files To/From the Cluster

### Linux/OS X
There are graphical tools to transfer files to and from the cluster but it is very handy to know the command line versions. For this we use the scp (secure copy) command.

The basic syntax is (and you need a terminal open on your local machine):

    :::bash
    scp username@server:remotefilename localfilename


The remote filename can be the full path or just the filename if it's in your home directory. The local filename can just be `.` for your current directory.

For example for a single file we can use:

    :::bash
    scp mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/pog.fa .


For a directory and its contents we can use the -r option :

    :::bash
    scp -r mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/ .

Common mistakes include:

* Forgetting the colon : separating the servername from the path
* Mistyping the path
* Forgetting the -r when copying directories
* Not including the destination filename/directory
* Running scp from the remote machine and not your local machine - check your command prompt people! Everyone does it at least once.

### Windows
The RC documentation has info on how to use filezilla here
[https://rc.fas.harvard.edu/resources/documentation/transferring-data/copying-data-to-and-fromodyssey-using-scp](https://rc.fas.harvard.edu/resources/documentation/transferring-data/copying-data-to-and-fromodyssey-using-scp)

<div style="color:red; margin: 50px">
Exercises
<ul>
<li>Transfer the /n/regal/informatics/workshops/Intermediate_Unix/Data/seq directory used above to your local machine. Try and get the command right first time and use pwd/cut and paste to minimize typing.</li>
</ul>
</div>


## 4. Redirecting Output and Backgrounding processes

### 4.1 Redirecting Output

We introduced redirecting output briefly in the Basic Unix workshop. We said that we could put the output from any command into a file using the > redirect operator. For instance:

    :::bash
    ls -l > mydir.dat

There are actually two types of output from a typical unix command stdout (standard out) which is what we were using above. There is also stderr (standard error) which is used for, not surprisingly, error messages. If we don't specify anything the shell assumes we're referring to stdout. If we want to differentiate them we refer to stdout as `1>` and stderr as `2>`. So if we want to put output and errors into separate files we'd do

    :::bash
    mycommand 1> mycmd.out 2> mycmd.err

If we have both types of output and want to put them into the same file (the most common action) we can do

    :::bash
    mycommand 2>&1 1> mycmd.out

This says to put the error output into the normal output stream and then the normal output stream into a file.

## 4.2 Backgrounding Processes

Most of the commands we've run in the previous workshop have finished almost immediately but with real analysis this is often not the case. To keep the command running but reclaim our command prompt we can end a command with & to tell the shell to run it in the background. 

For instance

    :::bash
    bash mylongcmd.sh > myfile.out &

If myfile.out exists this command will overwrite it. If you want to append to a file do:

    :::bash
    bash mylongcmd.sh >> myfile.out &


We can also push something into the background once it's running by using the ctrl-z command to suspend it and give us our command prompt back. Once it's suspended we can type bg to push it into the background and it'll start running again.

Other useful commands when you have background jobs are:

* `Ctrl-z bg` Suspend the job and then run it in the background
* `jobs` Lists all running jobs
* `kill %1` Kill job 1
* `kill %%` Kill all jobs
* `fg` Put the last job back into the foreground (or use fg n for the nth job)	

<div style="color:red; margin: 50px">
Exercises
<ul>

<li>Run the command 
  `bash /n/regal/informatics/workshops/Intermediate_Unix/Data/mylongscript.sh`
and put the output into a file</li>
<li>Suspend the job</li>
<li>Run the jobs command - what is it telling you?</li>
<li>Push the job into the background by typing bg</li>
<li>Rerun the jobs command - does the output make sense?</li>
<li>Kill your job</li>
</ul>
</div>


## 5. Slurm submission

We are now going to go over the basics of a slurm submission. You’re currently all logged into a login node where you can run short commands but anything that requires more than about 30 seconds of run time should be submitted to the slurm queue.

You can submit jobs completely on the command line but I recommend creating a script file and submitting that. This gives you a record of what you have run and the parameters you used. Here is a template for a typical slurm submit (sbatch) script.

    :::bash 
    #!/bin/bash
    #
    # These are comment lines starting with #
    #
    # These lines are interpreted by slurm.
    #SBATCH -J <jobname>        # The name of the job (can be any string – make is something readable)
    #SBATCH -N 1                # Ensure that all cores are on one machine
    #SBATCH -n <n>              # Use n cores for the job
    #SBATCH -t <n-nn:nn>        # Runtime in D-HH:MM
    #SBATCH -p <queuename>      # Partition to submit to (serial_requeue/general)
    #SBATCH --mem=<n>           # Memory pool for all cores in Mb (see also --mem-per-cpu)
    #SBATCH -o <outfile>.%A.out # File to which STDOUT will be written (%A is replaced by the jobid)
    #SBATCH -e <outfile>.%A.err # File to which STDERR will be written (%A is replaced by the jobid)
    #SBATCH --mail-type=<type>  # Type of email notification- BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=<myemail@harvard.edu> # Email to which notifications will be sent
    #SBATCH --constraint=holyib               # If you want to use the /n/holylfs storage
    
    # The command can use the command line  parameters $1 here
   
    source new-modules.sh
 
    module load <mymodule>
    
    some_command_here $1 > $1.out
    another_command < $1.out > $.new


In the sections below we'll go through the various pieces.

In practice it's handy to keep a template sbatch file handy so you can copy it when you want to run a new set of commands on the cluster.

### 5.1 The Header Line

    :::bash
    #!/bin/bash

This tells unix to use the /bin/bash command to execute the file.


### 5.2 The #SBATCH comment lines

In a bash script comment lines start with a #. The shell ignores these but when you submit a script to slurm all lines starting

    :::bash
    #SBATCH

are treated differently and tell slurm how to schedule your job.

### 5.3 The commands

Underneath the #SBATCH lines we can start running commands. If you’re using centrally installed software these can include module load commands.

### 5.4 Testing

Before submitting this script test this script on the command line by running it for a short while
using

    :::bash
    bash myscript.sh


If there is a typo or something wrong with the script it will fail almost immediately. Fix the problem and test until it runs.

### 5.5 Submitting
Submit your script to the cluster using

    :::bash
    sbatch myscript.h

When it’s submitted you’ll get a message containing the job id of the job.

### 5.6 Checking job status
If you’ve set the SBATCH parameters correctly you’ll get an email when your job runs/fails/etc. You can also check on the status of the job using the squeue command.

    :::bash
    squeue –u <username>


For an individual job you can see the status using

    :::bash
    sacct –j <jobid>


### 5.7 Killing jobs

If you have submitted something and you want to remove it from the queue use the scancel command

    :::bash
    scancel <jobid>

<div style="color:red; margin: 50px">
Exercises
<ul>
<li>Copy a template slurm submission script from
/n/regal/informatics/workshops/Intermediate_Unix/Slurm/template.sh to your home directory.</li>
<li>Edit this file and enter suitable values for all the parameters in the #SBATCH lines. We
recommend
<ul>
<li>Use 1 machine</li>
<li>Run on 1 core</li>
<li>Use 1Mbyte of memory</li>
<li>Allow for 5 minutes of run time</li>
<li>Submit to the serial_requeue queue</li>
</ul>
Inside the script run two commands `ls –l /tmp/` and `hostname` and append the output into a file.</li>
<li>Test your script on the command line (ctrl-c to get out of it)</li>
<li>Submit the script to slurm</li>
<li>Check the status of your job.</li>
</ul>
</div>


### Summary of SLURM commands
The table below shows a summary of SLURM commands, along with LSF equivalents and an example. These commands are described in more detail below along with links to the SLURM doc site.


<table style="margin: 20px; wrap=none">
<thead>
	<tr><th></th><th>SLURM COMMAND</th><th>SLURM EXAMPLE</th></tr>
</thead>
<tbody>
	<tr><td>Submit a batch serial job </td><td>sbatch</td><td>`sbatch runscript.sh`</td></tr>
	<tr><td>Run a script interatively</td><td>srun</td><td>`srun --pty -p interact -t 10 --mem 1000 /bin/hostname`</td></tr>
	<tr><td>Kill a job</td><td>scancel</td><td>`scancel 999999`</td></tr>
	<tr><td>View status of queues</td><td>squeue</td><td>`squeue -u akitzmiller`</td></tr>
	<tr><td>Check current job by id</td><td>sacct</td><td>`sacct -j 999999`</td></tr>
</tbody>
</table>


## 6. Searching for files – find

Find is a sophisticated recursive directory search which can locate files based on a pattern but also execute commands on files when you find them.

Basic syntax is :

    :::bash
    find <dir> -name <name pattern>


This will print all found files matching the pattern underneath the specified directory.

Example :

    :::bash
    find . -name “*.fastq”


The power comes when you add on the -exec option and append a command.

Example:

    :::bash
    find . -name “*.fastq” -exec ls -l {} \;


(You put {} where the filename should go and the command should always end with \;)

This command finds all fastq files and lists them.

Another usage for find is finding files that are newer or older than a certain time

    :::bash
    find . -mtime 7


Will find all files below the current directory modified more than 7 days ago. Use -7 for less than 7 days ago.

We can also use find for finding large files (or indeed small files)

    :::bash
    find . -size +1M 

And if we want to know the size of the files we can add an -exec option

    :::bash
    find . -size +1G -exec ls -ls {} \;

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find all files in your home directory older than 1 day</li>
<li> Find all files in /tmp/ newer than 7 days</li>
<li> Find all files under the /n/regal/informatics/workshops/Intermediate_Unix directory that end `.fa` and print their contents to the screen</li>
<li> Find all files under the /n/regal/informatics/workshops/Intermediate_Unix directory that are greater than 10G.  How big are they?</li>
</ul>
</div>

Extra: A good set of find examples is here [http://alvinalexander.com/unix/edu/examples/find.shtml](http://alvinalexander.com/unix/edu/examples/find.shtml)

## 6.1 Basic file searching - grep

The grep command looks for strings within files.  The basic use is :

   <pre>grep 'mystring' myfile</pre>

For example

   <pre>grep 'chr'  Data/AF1.bed </pre>

finds all lines in Data/AF1.bed that contain the string `chr`

(We'll come back to more advanced grep later)

## 7. Pipes – joining commands together using |

You can feed the output of one command to the input of another using the | character

For example you can run three commands one after the other 

    :::bash
    <cmd1> | <cmd2> | <cmd3>

This is best illustrated with some examples:

### Example 7.1

    :::bash
    ls –l | less

This is useful if you’re listing a big directory and the output won’t fit on one screen.

### Example 7.2

    :::bash
    find . -name “*.bed” | wc -l


This will find how many *.bed files there are below the current directory (wc –l returns the number of lines)

### Example 7.3

    :::bash
    find . -name “*.bed” -exec cat {} \; | wc -l


This finds how many lines total in all bed files found

### Example 7.4

    :::bash
    find . -name “*.bed” -exec cat {} \; |grep chr20 | wc -l


This finds how many chr20 lines there are in all bed files


### Example 7.5

Here's a more complicated but common and useful example

    :::bash
    cat myfile.dat |grep chr20 | awk ’{print $1}’ | sort | uniq –c > myfile.out


Here we have 5 different commands all chained together.   Breaking this down the command we have:


* lists the contents of myfile.dat (`cat`)
* for all lines that contain chr20  (`grep`)
* prints the string in the first column (`awk '{print $1}'`) 
* sorts the output (`sort`)
* prints only the unique strings and how many occurrences there are (`uniq -c`)
* puts the results into a file myfile.out (`>`)


Sometimes we have commands that absolutely need to take a file name rather than just read in input. In these cases we can use pipes and use the – character to replace the filename. For instance the samtools command on a file looks like :

    :::bash
    samtools view –b –S myfile.sam

If we want to feed input into samtools directly without using the myfile.sam file we need to do

    :::bash
    bowtie2 –x hg19 –U myfile.fq –p 32 | samtools view –v –S -

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find how many chr14 lines there are in
/n/regal/informatics/workshops/Intermediate_Unix/Data/AF2.bed and
/n/regal/informatics/workshops/Intermediate_Unix/Data/AF1.bed</li>
<li> Find and concatenate all .bed files under /n/regal/informatics/workshops/Intermediate_Unix/
and use pipe and another command to find the 10th row containing chr20 </li>
</ul>
</div>

## 8. The sort command

When we create data output files we often want to manipulate the contents by sorting. The unix `sort` command can be used to sort files very easily and in many different ways.

When we have columns of data we often want to sort on a column to find the highest or lowest entry. A typical command looks like:

    :::bash
    sort -nk4 <file> |less


There are a lot of useful options to sort.  The most useful are :


* `n` – sort numerically
* `k4` – sort starting on the 4th column
* `k4,5` – sort using the 4th and 5th columns only
* `r` – reverse sort
* `u` – sort and report unique lines
* `t”,”` – set the field delimiter to a comma



### Example 8.1

    :::bash
    sort -nk2 AF1.bed # sorts the file by the 2nd column


We can have multiple column options

    :::bash
    sort -k1,1 -k2,2n AF1.bed # sorts the file by chromosome first and then start coord


<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find the highest scoring 10 entries in the AF1.bed file (score is the 5th column)</li>
<li> Concatenate the AF1.bed and AF2.bed files and sort the results by chromosome and then score.</li>
<li> What is the 3rd lowest score in chr9 in the resulting file?</li>
</ul>
</div>

Extra:
More sort examples are at [http://www.theunixschool.com/2012/08/linux-sort-command-examples.html](http://www.theunixschool.com/2012/08/linux-sort-command-examples.html)

### Sorting big files

The sort command by default uses the /tmp/ directory to store intermediate files as it’s sorting. For very large files this can fill up the /tmp/ directory and your sort will fail. You can change where sort keeps its temporary files by using the –T option to specify a different directory. For instance if you want to use your current directory to store temp files use :

    :::bash
    sort –T . mybigfile.dat > mybigfile.sort


Of course you still have to have enough space in the new directory (use the `df –h` command to check)

## 9. Searching for strings – grep

`grep` is a fantastic command for searching through files and directories. The basic syntax is :

    :::bash
    grep <pattern> <file>


So to find all entries for chr20 in our AF1.bed file we’d do:

    :::bash
    grep ’chr20’ /n/regal/informatics/workshops/Intermediate_Unix/Data/AF1.bed


Useful options :

* `-v`   search for everything but the pattern
* `-n`   show the line number of the line found
* `-c`   show the count of the number of matched lines
* `-C`   n Show n lines context
* `-r`   search recursively down the directory tree (<file> is directory here)
* `-i`   ignore case
* `-H`   print the filename along with the file found
* `-f`   find files only
* `-d`   find directories only
* `-l`   only print the filename and not the line found (useful when there are multiple matches per file)


Extra :

Information about more complicated searching using regular expressions and egrep can be found here:
[http://ryanstutorials.net/linuxtutorial/grep.php](http://ryanstutorials.net/linuxtutorial/grep.php)

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find how many bed peaks (rows) there are in the AF2.bed file for chr10</li>
<li> Find which files contain the string AAAAAAAA in the Data directory</li>
</ul>
</div>


## 10. Manipulating file contents – awk

This is where things really get powerful. Awk is a ‘pattern scanning and processing’ utility. It lets you search and filter files based on column and pattern. The basic syntax is:

    :::bash
    awk ’pattern { action }’ filename


Again this is best shown by example

    :::bash
    awk ’$1 == “chr1” { print $1,$2,$3}’ AF1.bed


Here the pattern is `$1 == “chr1”` or column1 = chr1 and the action is print columns 1,2 and 3

We can do more interesting things with the pattern e.g.

    :::bash
    awk ’$2 > 1000000 && $3 < 2000000 { print $0}’ AF1.bed


This only prints lines where the region is between 1 and 2Mb. The `$0` represents the whole line

We could also do this by missing out the whole action

    :::bash
    awk ’$2 > 1000000 && $3 < 2000000’ AF1.bed


We can also search for substrings using the ~ character

    :::bash
    awk ’$1 ~ /1/ { print $1}’ AF1.bed


only prints out lines where the first column contains a 1

Similarly

    :::bash
    awk ’$1 !~ /1/ {print $1}’

This prints the first column where the first column doesn't contain a 1

There are other things you can include

    :::bash
    awk ’$1 ~ /^1/ {print $2}’


This prints the 2nd column where the 1st column starts with a 1
You can reference the line number and number of fields using NR and NF.

    :::bash
    awk ’NR > 1’ # will only print out rows 2 - end
    
    awk ’NF == 12’ # will only print out rows with exactly 12 fields
     
    awk ’NR % 4 == 1’ # will only print lines where the line_number / 4 has a remainder of 1


### Combining awk with sort and uniq

The uniq command omits repeated lines. It is often used to count repeated lines using the -c option

For example:

    :::bash
    awk ’{print $1}’ Data/AF1.bed |sort | uniq -c


This counts how often each chromosome appears in the bed file.

There’s much more to awk but these commands will get you a long way. Let’s now do something useful :

<div style="color:red; margin: 50px">
Exercises:
<ul>
<li>Find how many entries there are in the AF1.bed file for each chromosome</li>
<li>Find how many entries scoring > 300 there are in the AF1.bed file for each chromosome</li>
<li>Find the read lengths in the fastq file (Hint: each fastq entry has 4 lines and the read length is on the 1st line of the entry)</li>
</ul>
</div>

