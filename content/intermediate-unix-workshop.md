Title: Intermediate Unix Workshop
Date: 2016-01-13 16:05
Category: Tutorials
Author: mclamp@harvard.edu
Summary: This tutorial is intended for people who are familiar with the basics of unix but want to learn more about manipulating files and running commands.


## Notes

This text is available at<br><Br>


*   [http://informatics.fas.harvard.edu/intermediate-unix-workshop](http://informatics.fas.harvard.edu/intermediate-unix-workshop/)


## 1. Introduction and prerequisites

This tutorial is intended for people who are familiar with the basics of unix but want to learn more
about manipulating files and running commands.

*   A Harvard FAS RC cluster account 
    [https://account.rc.fas.harvard.edu/request/](https://account.rc.fas.harvard.edu/request/)

*   A terminal program so you can log into the cluster 
    [https://rc.fas.harvard.edu/resources/accessand-login/](https://rc.fas.harvard.edu/resources/accessand-login/)

*   The ability to log into the cluster 
    [https://rc.fas.harvard.edu/resources/access-and-login/](https://rc.fas.harvard.edu/resources/access-and-login/)

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
* perl regex
* bash for loops

### Functionality Covered

We call this the intermediate tutorial and it is intended for people who are comfortable with creating
files and moving around the filesystem but want to move onto more manipulation of data files.

## 3.  Transferring Files To/From the Cluster

### Linux/OS X
There are graphical tools to transfer files to and from the cluster but it is very handy to know the
command line versions. For this we use the scp (secure copy) command.

The basic syntax is (and you need a terminal open on your local machine):
<pre style="margin-top: 20px">
scp username@server:remotefilename localfilename
</pre>

The remote filename can be the full path or just the filename if it's in your home directory. The local
filename can just be . for your current directory.

For example for a single file we can use:

<pre style="margin-top: 20px">
scp mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/pog.fa .
</pre>

For a directory and its contents we can use the -r option :

<pre style="margin-top: 20px">
scp -r mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/ .
</pre>

Common mistakes include:

<br>
<ul>
<li>Forgetting the colon : separating the servername from the path</li>
<li>Mistyping the path</li>
<li>Forgetting the -r when copying directories</li>
<li>Not including the destination filename/directory</li>
<li>Running scp from the remote machine and not your local machine - check your command<
prompt people! Everyone does it at least once.</li>
</ul>

### Windows
The RC documentation has info on how to use winscp here
[https://rc.fas.harvard.edu/resources/documentation/transferring-data/copying-data-to-and-fromodyssey-using-scp](https://rc.fas.harvard.edu/resources/documentation/transferring-data/copying-data-to-and-fromodyssey-using-scp)

<div style="color:red; margin: 50px">
Exercises
<ul>
<li>Transfer the /n/regal/informatics/workshops/Intermediate_Unix/Data/seq directory used above
to your local machine. Try and get the command right first time and use pwd/cut and paste to
minimize typing.</li>
</ul>
</div>


## 4. Redirecting Output and Backgrounding processes

###4.1 Redirecting Output

We introduced redirecting output briefly in the Basic Unix workshop. We said that we could put the
output from any command into a file using the > redirect operator. For instance:

<pre style="margin-top: 20px">
ls -l > mydir.dat
</pre>

There are actually two types of output from a typical unix command stdout (standard out) which is
what we were using above. There is also stderr (standard error) which is used for, not surprisingly,
error messages. If we don't specify anything the shell assumes we're referring to stdout. If we want to
differentiate them we refer to stdout as 1> and stderr as 2>. So if we want to put output and errors
into separate files we'd do

<pre style="margin-top: 20px">
mycommand 1> mycmd.out 2> mycmd.err
</pre>

If we have both types of output and want to put them into the same file (the most common action) we
can do

<pre style="margin-top: 20px">
mycommand 2>&1 1> mycmd.out
</pre>
This says to put the error output into the normal output stream and then the normal output stream
into a file.

## 4.2 Backgrounding Processes

Most of the commands we've run in the previous workshop have finished almost immediately but with real analysis this
is often not the case. To keep the command running but reclaim our command prompt we can end a
command with & to tell the shell to run it in the background. 

For instance

<pre style="margin-top: 20px">
bash mylongcmd.sh > myfile.out &
</pre>

If myfile.out exists this command will overwrite it. If you want to append to a file do:

<pre style="margin-top: 20px">
bash mylongcmd.sh >> myfile.out &
</pre>

We can also push something into the background once it's running by using the ctrl-z command to
suspend it and give us our command prompt back. Once it's suspended we can type bg to push it into
the background and it'll start running again.

Other useful commands when you have background jobs are:


<ul>
<li>Ctrl-z bg Suspend the job and then run it in the background</li>
<li>jobs Lists all running jobs</li>
<li>kill %1 Kill job 1</li>
<li>kill %% Kill all jobs</li>
<li>fg Put the last job back into the foreground (or use fg n for the nth job)</li>
</ul>

<div style="color:red; margin: 50px">
Exercises
<ul>

<li>Run the command 
  <pre>bash /n/regal/informatics/workshops/Intermediate_Unix/Data/mylongscript.sh</pre> and put the output
into a file</li>
<li>Suspend the job</li>
<li>Run the jobs command - what is it telling you?</li>
<li>Push the job into the background by typing bg</li>
<li>Rerun the jobs command - does the output make sense?</li>
<li>Kill your job</li>
</ul>
</div>


## 5. Slurm submission

We are now going to go over the basics of a slurm submission. You’re currently all logged into a login
node where you can run short commands but anything that requires more than about 30 seconds of
run time should be submitted to the slurm queue.

You can submit jobs completely on the command line but I recommend creating a script file and
submitting that. This gives you a record of what you have run and the parameters you used. Here is
a template for a typical slurm submit (sbatch) script.

<pre style="margin-top: 20px">
#!/bin/bash
#
# These are comment lines starting with #
#
# These lines are interpreted by slurm.
#SBATCH -J <jobname> # The name of the job (can be any string – make is something readable)
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n <n> # Use n cores for the job
#SBATCH -t <n-nn:nn> # Runtime in D-HH:MM
#SBATCH -p <queuename> # Partition to submit to
#SBATCH --mem=<n> # Memory pool for all cores in Mb (see also --mem-per-cpu)
#SBATCH -o <outfile>.%A.out # File to which STDOUT will be written (%A is replaced by the jobid)
#SBATCH -e <outfile>.%A.err # File to which STDERR will be written (%A is replaced by the jobid)
#SBATCH --mail-type=<type> # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<myemail@harvard.edu> # Email to which notifications will be sent

# The command can use the input parameter $1 here

module load <mymodule>

some_command_here $1 > $1.out
another_command < $1.out > $.new


</pre>

In the sections below we'll go through the various pieces.

### 5.1 The Header Line
<pre style="margin-top: 20px">
\#!/bin/bash
</pre>
This tells unix to use the /bin/bash command to execute the file.


### 5.2 The #SBATCH comment lines

In a bash script comment lines start with a #. The shell ignores these but when you submit a script
to slurm all lines starting
<pre style="margin-top: 20px">
\#SBATCH
</pre>
are treated differently and tell slurm how to schedule your job.

### 5.3 The commands

Underneath the #SBATCH lines we can start running commands. If you’re using centrally installed
software these can include module load commands.

### 5.4 Testing

Before submitting this script test this script on the command line by running it for a short while
using
<pre style="margin-top: 20px">
bash myscript.sh
</pre>

If there is a typo or something wrong with the script it will fail almost immediately. Fix the problem
and test until it runs.

### 5.5 Submitting
Submit your script to the cluster using
<pre style="margin-top: 20px">
sbatch myscript.h
</pre>
When it’s submitted you’ll get a message containing the job id of the job.

### 5.6 Checking job status
If you’ve set the SBATCH parameters correctly you’ll get an email when your job runs/fails/etc. You
can also check on the status of the job using the squeue command.

<pre style="margin-top: 20px">
squeue –u username
</pre>

For an individual job you can see the status using
<pre style="margin-top: 20px">
sacct –j jobid
</pre>

### 5.7 Killing jobs

If you have submitted something and you want to remove it from the queue use the scancel command
<pre style="margin-top: 20px">
scancel <jobid>
</pre>

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
Run two commands ls –l /tmp/ and hostname and append the output into a file.</li>
<li>Test your script on the command line (ctrl-c to get out of it)</li>
<li>Submit the script to slurm</li>
<li>Check the status of your job.</li>
</ul>
</div>


### Summary of SLURM commands
The table below shows a summary of SLURM commands, along with LSF equivalents and an
example. These commands are described in more detail below along with links to the SLURM doc
site.


<table style="border: 1px solid black; margin: 20px; wrap=none"><tbody><th><td>SLURM COMMAND</td><td>SLURM EXAMPLE</td></th>
<tr><td>Submit a batch serial job </td><td>sbatch</td><td>sbatch runscript.sh</td></tr>
<tr><td>Run a script interatively</td><td>srun</td><td>srun --pty -p interact -t 10 --mem 1000 /bin/hostname</td></tr>
<tr><td>Kill a job</td><td>scancel</td><td>scancel 999999</td></tr>
<tr><td>View status of queues</td><td>squeue</td><td>squeue -u akitzmiller</td></tr>
<tr><td>Check current job by id</td><td>sacct</td><td>sacct -j 999999</td></tr>
</tbody></table>


## 6. Searching for files – find

Find is a sophisticated recursive directory search which can locate files based on a pattern but also
execute commands on files when you find them.

Basic syntax is :

<pre style="margin-top: 20px">
find \<di\r> -name \<name pattern\>
</pre>

This will list all found files.

Example :

<pre style="margin-top: 20px">
find . -name “*.fastq”
</pre>

The power comes when you add on the -exec option and append a command.

Example:

<pre style="margin-top: 20px">
find . -name “*.fastq” -exec ls -l {} \;
</pre>

(You put {} where the filename should go and the command should always end with \;)

This command files all fastq files and lists them.

Another usage for find is finding files that are newer or older than a certain time

<pre style="margin-top: 20px">
find . -mtime 7
</pre>

Will find all files below the current directory modified more than 7 days ago. Use -7 for less than 7 days ago.

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find all files in your home directory older than 1 day</li>
<li> Find all files in /tmp/ newer than 7 days</li>
<li> Find all files under the /n/regal/informatics/workshops/Basic_UNIX directory that contain the string chr14</li>
</ul>
</div>

Extra: A good set of find examples is here [http://alvinalexander.com/unix/edu/examples/find.shtml](http://alvinalexander.com/unix/edu/examples/find.shtml)


## 7. Pipes – joining commands together using |

You can feed the output of one command to the input of another using the | character

For example you can run three commands one after the other 

<pre style="margin-top: 20px">
\<cmd1\> | \<cmd2\> | \<cmd3\>
</pre>

This is best illustrated with some examples:

### Example 7.1
<pre style="margin-top: 20px">
ls –l | less
</pre>

This is useful if you’re listing a big directory and the output won’t fit on one screen.

### Example 7.2
<pre style="margin-top: 20px">
find . -name “*.bed” | wc -l
</pre>

This will find how many *.bed files there are below the current directory (wc –l returns the number of lines)

### Example 7.3
<pre style="margin-top: 20px">
find . -name “*.bed” -exec cat {} \; | wc -l
</pre>

This finds how many lines total in all bed files found

### Example 7.4
<pre style="margin-top: 20px">
find . -name “*.bed” -exec cat {} \; |grep chr20 | wc -l
</pre>

This finds how many chr20 lines there are in all bed files


### Example 7.5

Here's a more complicated but common and useful example

<pre style="margin-top: 20px">
cat myfile.dat |grep chr20 | awk ’{print $1}’ | sort | uniq –c > myfile.out
</pre>

Here we have 5 different commands all chained together.   Breaking this down the command we have:



<ul>
<li> lists the contents of myfile.dat (cat)</li>
<li> for all lines that contain chr20  (grep)</li>
<li> prints the string in the first column (awk '{print $1}') </li>
<li> sorts the output (sort)</li>
<li> prints only the unique strings and how many occurrences there are (uniq -c)</li>
<li> puts the results into a file myfile.out (>)</li>
</ul>


Sometimes we have commands that absolutely need to take a file name rather than just read in input.
In these cases we can use pipes and use the – character to replace the filename. For instance the
samtools command on a file looks like :

<pre style="margin-top: 20px">
samtools view –b –S myfile.sam
</pre>
If we want to feed input into samtools directly without using the myfile.sam file we need to do

<pre style="margin-top: 20px">
bowtie2 –x hg19 –U myfile.fq –p 32 | samtools view –v –S -
</pre>

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find how many chr14 lines there are in
/n/regal/informatics/workshops/Intermediate_Unix /Data/AF2.bed and
/n/regal/informatics/workshops/Intermediate_Unix /Data/AF1.bed</li>
<li> Find and concatenate all .bed files under /n/regal/informatics/workshops/Intermediate_Unix/
and use pipe and another command to find the 10th row containing chr20 </li>
</ul>
</div>

## 7. The sort command

When we create data output files we often want to manipulate the contents by sorting. The unix
sort command can be used to sort files very easily and in many different ways.

When we have columns of data we often want to sort on a column to find the highest or
lowest entry. A typical command looks like:

<pre style="margin-top: 20px">
sort -nk4 <file> |less
</pre>

There are a lot of useful options to sort.  The most useful are :


<ul>
<li>n – sort numerically</li>
<li>k4 – sort starting on the 4th column</li>
<li>k4,5 – sort using the 4th and 5th columns only</li>
<li>r – reverse sort</li>
<li>u – sort and report unique lines</li>
<li>t”,” – set the field delimiter to ,</li>
</ul>


### Example 7.1
<pre style="margin-top: 20px">
sort -nk2 AF1.bed # sorts the file by the 2nd column
</pre>

We can have multiple column options

<pre style="margin-top: 20px">
sort -k1,1 -k2,2n AF1.bed # sorts the file by chromosome first and then start coord
</pre>

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find the highest scoring 10 entries in the AF1.bed file (score is the 5th column)</li>
<li> Concatenate the AF1.bed and AF2.bed files and sort the results by chromosome and then score.</li>
<li> What is the 3rd lowest score in chr9 in the resulting file?</li>
</ul>
</div>

Extra:
More sort examples are at http://www.theunixschool.com/2012/08/linux-sort-command-examples.html

### Sorting big files

The sort command by default uses the /tmp/ directory to store intermediate files as it’s sorting. For
very large files this can fill up the /tmp/ directory and your sort will fail. You can change where sort
keeps its temporary files by using the –T option to specify a different directory. For instance if you
want to use your current directory to store temp files use :

<pre style="margin-top: 20px">
sort –T . mybigfile.dat > mybigfile.sort
</pre>

Of course you still have to have enough space in the new directory (use the df –h command to check)

## 8. Searching for strings – grep

Grep is a fantastic command for searching through files and directories. The basic syntax is :

<pre style="margin-top: 20px">
grep \<pattern> \<file>
</pre>

So to find all entries for chr20 in our AF1.bed file we’d do:

<pre style="margin-top: 20px">
grep ’chr20’ /n/regal/informatics/workshops/Intermediate_UNIX/Data/AF1.bed
</pre>

Useful options :
<ul>
<li> -v   search for everything but the pattern</li>
<li> -n   show the line number of the line found</li>
<li> -c   show the count of the number of matched lines</li>
<li> -C   n Show n lines context</li>
<li> -r   search recursively down the directory tree (<file> is directory here)</li>
<li> -i   ignore case</li>
<li> -H   print the filename along with the file found</li>
<li> -f   find files only</li>
<li> -d   find directories only</li>
<li> -l   only print the filename and not the line found (useful when there are multiple matches per file)</li>
</ul>

Extra :

Information about more complicated searching using regular expressions and egrep can be found here:
http://ryanstutorials.net/linuxtutorial/grep.php

<div style="color:red; margin: 50px">
Exercises :
<ul>
<li> Find how many bed peaks (rows) there are in the AF2.bed file for chr10</li>
<li> Find which files contain the string AAAAAAAA in the Data directory</li>
</ul>
</div>


## 9. Manipulating file contents – awk

This is where things really get powerful. Awk is a ‘pattern scanning and processing’ utility. It lets
you search and filter files based on column and pattern. The basic syntax is:

<pre style="margin-top: 20px">
awk ’pattern { action }’ filename
</pre>

Again this is best shown by example

<pre style="margin-top: 20px">
awk ’$1 == “chr1” { print $1,$2,$3}’ AF1.bed
</pre>

Here the pattern is $1 == “chr1” or column1 = chr1 and the action is print columns 1,2 and 3

We can do more interesting things with the pattern e.g.

<pre style="margin-top: 20px">
awk ’$2 > 1000000 && $3 < 2000000 { print $0}’ AF1.bed
</pre>

This only prints lines where the region is between 1 and 2Mb. The $0 represents the whole line

We could also do this by missing out the whole action

<pre style="margin-top: 20px">
awk ’$2 > 1000000 && $3 < 2000000’ AF1.bed
</pre>

We can also search for substrings using the ~ character

<pre style="margin-top: 20px">
awk ’$1 ~ /1/ { print $1}’ AF1.bed
</pre>

only prints out lines where the first column contains a 1

Similarly

<pre style="margin-top: 20px">
awk ’$1 !~ /1/ {print $1}’
</pre>

This prints the first column where the first column doesnt’ contain a 1

There are other things you can include

<pre style="margin-top: 20px">
awk ’$1 ~ /^1/ {print $2}’
</pre>

This prints the 2nd column where the 1st column starts with a 1
You can reference the line number and number of fields using NR and NF.

<pre style="margin-top: 20px">
awk ’NR > 1’ # will only print out rows 2 - end
</pre>
<pre style="margin-top: 20px">
awk ’NF == 12’ # will only print out rows with exactly 12 fields
</pre>
<pre style="margin-top: 20px">
awk ’NR % 4 == 1’ # will only print lines where the line_number / 4 has a remainder of 1
</pre>

### Combining awk with sort and uniq

The uniq command omits repeated lines. It is often used to count repeated lines using the -c option

For example:
<pre style="margin-top: 20px">
awk ’{print $1}’ Data/AF1.bed |sort | uniq -c
</pre>

This counts how often each chromosome appears in the bed file.

There’s much more to awk but these commands will get you a long way. Let’s now do something
useful :
<div style="color:red; margin: 50px">
Exercises:
<ul>
<li>Find how many entries there are in the AF1.bed file for each chromosome</li>
<li>Find how many entries scoring > 300 there are in the AF1.bed file for each chromosome</li>
<li>Find the read lengths in the fastq file (Hint: each fastq entry has 4 lines and the read length is on the 1st line of the entry)</li>
</ul>
</div>


## 10. Perl regular expressions

Perl is a complete programming language but there is one aspect of it that is really useful on the command
line. This is it’s ability to search and replace strings on the fly.

An example of the syntax is:

<pre style="margin-top: 20px">
perl -pe ’s/sausage/melon/g’ myfile
</pre>

This will take each line of myfile and change each instance (the final g) of sausage into melon. 

More powerfully we can have wildcards:

<pre style="margin-top: 20px">
perl -pe ’s/.*Length=//’ myfile
</pre>

Let's break the search pattern .*Length= down :

*  The . means match any character and the * means match any number of characters.
*  Length= means find Length= in each line

Now let's look at the replace pattern.  In this case it's empty so what we have is

 - match any characters up to and including Length= and replace them with nothing (i.e. delete them)


Let’s do something more complicated

<pre style="margin-top: 20px">
perl -pe ’s/.*:(.*?)/$1/’ myfile
</pre>

Again let's break this down.

The search pattern is .*:(.*?) which looks nasty but let's go through it

*   .*   we know means any number of characters
*   :    match a colon
*   (    the brackets mean remember this match for later
*   .*?  means match any number of characters but as small a number as possible to fit in with the other search elements

So we're searching and saving all characters after the last colon in the line

The replace pattern is just $1 which means replace with whatever was matched within the first set of brackets in the search pattern.

So this removes everything before the final colon.

Matches can be specified more precisely than . and.*

We can use:


</ul>
<li>\S+ for one or more non-whitespace (* is 0 or more and + is 1 or more)</li>
<li>\s+ for whitespace</li>
<li>\d+ for digits</li>
<li>\w+ for words</li>
<li>^ Start of line</li>
<li>$ end of line</li>
<li>\n newline char</li>
<li>\t tab char</li>
<li>\. Will match a . character (in general \ will escape a character)</li>
</ul>

This will reverse the first and 2nd words
<pre style="margin-top: 20px">
perl -pe ’s/^(\w+)\s+(\w+)/$2 $1/’ myfile
</pre>

And of course we can use pipes

<pre style="margin-top: 20px">
awk ’$5 > 200 { print $1}’ |perl -pe ’s/chr//’
</pre>

This will print out all the column 1 chromosome labels for scores > 200 and strip of the chr string

<div style="color:red; margin: 50px">
Exercises:
<ul>
<li>1. Convert a fastq file to fasta format using one line. Use the
/n/regal/informatics/workshops/Intermediate_Unix/Data/SRR866428_1.fastq file.
Note:
Fastq format looks like
<pre style="margin-top: 20px">
Header line
GGTTATTAGGGTGGCAGAGCCAGGAAATTGCGT
+some stuff here
Quality line
</pre>

Fasta format looks like
<pre style="margin-top: 20px">
>header line here
GGTTATTAGGGTGGCAGAGCCAGGAAATTGCGT
>another header line
CCTCTAAGGCGGGCCACTGTGCCAAATTCTCTA
</li>
<li>2. Extract the index strings from the fastq file and estimate the frequency distribution
(Use the /n/regal/informatics/workshops/Intermediate_Unix/DF_2.R1.fastq file)
Index strings are on the end of the header line and look like CGTACTAG (or a similar length sequence)</li>
</ul>
</div>

Extra: Useful information about perl in command lines
http://www.softpanorama.org/Scripting/Perlorama/perl_in_command_line.shtml

