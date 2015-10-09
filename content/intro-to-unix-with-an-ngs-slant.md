Title: Intro to Unix with an NGS slant
Date: 2015-01-01 11:00
Category: Tutorials

This tutorial introduces basic Linux operation in the context of many of the operations that might be done for NGS data  processing.

###  Notes

This text is available at
*   [http://informatics.fas.harvard.edu/intro-to-unix-with-an-ngs-slant](http://informatics.fas.harvard.edu/intro-to-unix-with-an-ngs-slant/)


The input files can be found at
*   /n/regal/informatics/workshops/Basic_Unix/Data



The sbatch files will appear on the filesystem at
*   /n/regal/informatics/workshops/Basic_Unix/Slurm



The output files will appear on the filesystem at
*   /n/regal/informatics/workshops/Basic_Unix/Output



### Prerequisites

*   A Harvard FAS RC cluster account
*   The ability to log into the cluster using the command line

Both of these are covered in the first 14 slides of the RC Intro to Unix slides here

[https://software.rc.fas.harvard.edu/training/intro_unix/latest/#(1)](https://software.rc.fas.harvard.edu/training/intro_unix/latest/#(1))

This is a very good intro put together by John Brunelle and is worth a read.  We’ll be covering most of the things but with more of an emphasis on file content manipulation.



### **Logging in**

From an OS X Terminal  (or an xterm or similar) type





<pre>ssh -X [myusername@login.rc.fas.harvard.edu](mailto:myusername@login.rc.fas.harvard.edu)</pre>

You’ll be prompted for your password and then for your openauth (JAuth) token.  If all is successful then you’ll see something like

<pre>[mclamp@Micheles-iMac:~]$ ssh [mclamp@login.rc.fas.harvard.edu](mailto:mclamp@login.rc.fas.harvard.edu)
Odyssey2
Login issues? See [https://rc.fas.harvard.edu/resources/support/](https://rc.fas.harvard.edu/resources/support/)

Password:
Verification code:
Last login: Tue Apr 28 12:09:14 2015 from 10.255.12.32

[mclamp@rclogin03 ~]$

</pre>

You are now logged into the cluster and can start typing in commands.

(Aside:   If you’re used to solely working on your laptop or desktop it’s easy to forget which terminal is on the cluster and which is on your laptop.  No good trying to submit jobs to the cluster from your laptop!   As you get used to command line work you’ll get into the habit of just looking at the command prompt before you type anything in

<pre>[mclamp@rclogin03 ~]$
</pre>

Usually this will remind you who you are and which machine you’re just about to run something on.

<span style="color: #ff0000;">Hands on:</span>

<span style="color: #ff0000;">Log into the cluster now from your machine.    Note the difference between the command prompt I get and the one you have.</span>

## The Basics

Although we’re going to go through the basics these few commands will form most of the activity you have on the cluster.  They’re all to do with files/directories/copying/deleting etc.  Let’s start.

## 1.  Directories

### 1.1  Paths

Let’s dive in with some commands right away.    When you log in you land in your home directory.  You can find out where this is by typing pwd

<pre>mclamp@rclogin03 ~]$ pwd
/n/home_rc/mclamp
[mclamp@rclogin03 ~]$

</pre>

On unix directories are separated by a forward slash (unlike windows) .   Everything hangs off the root directory [ / ] so regardless of which physical disk or filesystem your files live on full paths always start with [ / ]





General Note: Typing man <command> will bring up the help page for any command.  This will list what it does, definitions of the different options and, if you’re lucky, some examples.

<span style="color: #ff0000;">Hands on:</span>



<span style="color: #ff0000;">1  Try the pwd command in your home directory now.   Note the differences.</span>



<span style="color: #ff0000;"> </span>





### <span style="color: #000000;">1.2  Changing directories</span>

The cd command changes directories.  You can either enter the full path

<pre>cd /n/regal/informatics/workshops
</pre>

or you can navigate relative to your current directory with .. which means move up a directory.   The equivalent to the previous command using this notation would be

<pre>cd ../../regal/informatics/workshops
</pre>





Note:   Just typing cd with no argument will take you back to your home directory.   You can also reference your home directory using the ~ character.

### 1.3 Listing directories and files



Knowing what’s in directories is useful and to do this we use the ls command.  ls takes many options; by far the most used by me is

<pre>ls -ltra

</pre>

*   l – long form
*   t – sort by time
*   r – reverse list
*   a – list all files

If you type this in a directory you’ll get a full listing of all the files with the most recent ones at the bottom.  This makes it easy to see the most recently changed files in one command.

You can use all, none or any combination of these options.  Other useful ones include -F   -s  and -h

The  output of ls which includes -l will give the long form of a directory listing which looks like

<pre>[mclamp@rclogin03 Mouse_ATAC]$ ls -ltra
total 36875604
drwxrwxrwx 4 mclamp rc_admin 4096 Jun 17 17:58 ..
-rw-r--r-- 1 mclamp rc_admin 4163383591 Jun 17 18:20 Ctrl1_atac_SP2.R1.fastq
-rw-r--r-- 1 mclamp rc_admin 4163383591 Jun 17 18:21 Ctrl1_atac_SP2.R2.fastq
-rw-r--r-- 1 mclamp rc_admin [4151176918](http://nrs.harvard.edu/urn-3:hul.ois:primoDeepSearch?institution=HVD&vid=HVD&tab=everything&search_scope=everything&mode=Basic&onCampus=false&displayMode=full&highlight=true&pcAvailabiltyMode=true&bulkSize=30&query=isbn%2Ccontains%2C4151176918&displayField=all) Jun 17 18:23 Ctrl2_atac_SP2.R1.fastq
-rw-r--r-- 1 mclamp rc_admin [4151176918](http://nrs.harvard.edu/urn-3:hul.ois:primoDeepSearch?institution=HVD&vid=HVD&tab=everything&search_scope=everything&mode=Basic&onCampus=false&displayMode=full&highlight=true&pcAvailabiltyMode=true&bulkSize=30&query=isbn%2Ccontains%2C4151176918&displayField=all) Jun 17 18:23 Ctrl2_atac_SP2.R2.fastq
-rw-r--r-- 1 mclamp rc_admin 5531764835 Jun 17 18:25 LPS1_atac_SP2.R1.fastq
-rw-r--r-- 1 mclamp rc_admin 5531764835 Jun 17 18:26 LPS1_atac_SP2.R2.fastq
drwxr-xr-x 2 mclamp rc_admin 4096 Jun 17 18:26 .
-rw-r--r-- 1 mclamp rc_admin 5033941307 Jun 17 18:27 LPS2_atac_SP2.R1.fastq
-rw-r--r-- 1 mclamp rc_admin 5033941307 Jun 17 18:28 LPS2_atac_SP2.R2.fastq

The columns are permissions, ?,owner, group, size in bytes, date modified, filename

The . file means the current directory (as opposed to .. which means the parent directory)

The permissions look complicated but are mostly straightforward.  

**d**rwxrwxrwxrwx     A d at the front means it's a directory

-**rw-**r--r--       The next three are user permissions (this case read and write)

-rw-**r--**r--       The second three are group permissions (read)

-rw-r--**r--** The final three are world permissions (read)

Directories have to have read and execute permissions for people to read them and read,write and execute for people to write.  Hence drwxrwxrwx means everybody and read and write to this directory.
</pre>

Note:  Command line editing

*   **ctrl-a   go to beginning of line**
*   **ctrl-e   got to end of line**
*   **tab tab  will give a list of all available files or tab completing (best used by starting typing and then tab tab**
*   **ctrl-f   forward a char**
*   **ctrl-b  back a char**
*   **ctrl-d   delete a char**
*   **ctrl-k   delete to end of line (ctrl-a ctrl-k is popular)**
*   ctrl-u   delete to beginning of line
*   esc-f    forward a word
*   esc-b   back a word
*   esc-d   delete a word

More useful shortcuts

*   ctrl-c   cancel current command
*   ctrl-s  stops output to screen
*   ctrl-p  restarts output to screen
*   ctrl-l    clears screen
*   ctrl-z   put current command in the background
*   fg        puts most recent backgorund command in foreground
*   jobs    lists backgrounded jobs  (kill %1 kills first backgrounded job,  kill %% kills all jobs)

More shortcuts here [http://www.skorks.com/2009/09/bash-shortcuts-for-maximum-productivity/](http://www.skorks.com/2009/09/bash-shortcuts-for-maximum-productivity/)

<span style="color: #ff0000;">Hands on :</span>

1.  <span style="color: #ff0000;">cd into the workshop data directory /n/regal/informatics/workshops/Basic_Unix</span>
2.  <span style="color: #ff0000;">Using a combination of cd and ls with different options investigate the contents of this directory.   Which options do you find most useful?</span>
3.  <span style="color: #ff0000;">The df command lists all the different filesystems, where they are mounted, how big they are and how full they are.   Try this now and use the man page to work out what the columns mean.</span>
4.  <span style="color: #ff0000;">Play with the command line editing commands until you can navigate a command line quickly</span>

<span style="color: #ff0000;"> </span>



### 1.4 Making and removing directories

The mkdir and rmdir make and remove directories.   Note that rmdir can only be used if the directory is empty.





<pre>mkdir pogdir
cd pogdir
ls -l 
cd ../
rmdir pogdir</pre>

You often want to remove a directory and its contents use rm -rf for this

<pre>mkdir pogdir
cd pogdir
ls -l 
cd ../
rm -rf pogdir</pre>



Be very careful of rm -rf !!

<span style="color: #ff0000;">Hands on:</span>

<span style="color: #ff0000;">1\. Make a user directory for yourself (use your username) under</span>  <span style="color: #ff0000;">/n/regal/informatics/workshops/Basic_Unix/Users</span>

<span style="color: #ff0000;"> </span>





## 2.  Files


### 2.1 Moving, copying and deleting files


#### Moving files

<pre>mv myfile mynewfile</pre>

You can use the ../ and

<pre>mv myfile ../
my *.fastq ../</pre>

<pre>mv myfile ../mynewfile</pre>

(also works with directories)

#### Copying files

<pre>cp myfile mynewfile</pre>

<pre>cp myfile ../
cp *.bam ../</pre>

<pre>cp myfile ../mynewfile</pre>

#### Removing files

<pre>rm myfile
rm *.bai

</pre>

rm -i myfile  will give you a prompt for each file

#### Copying directories

Usually you want to copy the directory and its contents.  You need to add the -r option to copy recursively

<pre>cp -r mydir mynewdir</pre>





If you want to copy the directory and keep all the timestamps the same your can use the -P option with cp

<pre>cp -P mydir mynewdir</pre>

or <span style="font-size: 1rem; line-height: 1;">rsync</span>





<pre>rsync -avz mydir mynewdir</pre>



#### Linking Files



For large files you can link rather than copy (saves disk space).  The syntax is

<pre>ln -s <real file> <linked file></pre>





<span style="color: #ff0000;">Hands on:</span>



1.  <span style="color: #ff0000;">Copy the /n/regal/informatics/workshops/Basic_Unix/Data directory and its contents to your user directory.</span>
2.  <span style="color: #ff0000;">Find the largest file in your data directory,  remove it and make a symlink instead.</span>
3.  <span style="color: #ff0000;">Remove your Data directory and its contents</span>



### 2.2  Creating and looking at files





#### **Nano is a basic text editor**

A user-friendly editor to create files is nano.     When you type nano to start editing the window helpfully lists the available comands at the bottom of the screen so you don’t have to remember anything.



#### **Other methods**

<pre>cat  filename  >  newfilename</pre>

<pre>echo  "some string here"   >  myfile</pre>

Sometimes you just want to have a peek at the start or end of a file.  The commands head and tail are very useful here

<pre>head myfile                    # will list the first 10 lines of a file</pre>

<pre>tail myfile                    # will list the last 10 lines of a file</pre>

If you want more lines

<pre>head -100 myfile               # The first 100 lines</pre>

<pre>tail -100 myfile               # The last 100 lines</pre>

One *very* useful command is tail -f <filename>



This lists the last 10 lines of a file and if the file gets bigger will also list any new lines.   Very useful if you’re running a command and want to see how the output file is growing.



<span style="color: #ff0000;">Hands on:</span>

<span style="color: #ff0000;">We have some sample data files in /n/regal/informatics/workshops/Basic_UNIX/Dat</span>

1.  <span style="color: #ff0000;">List that directory</span>
2.  <span style="color: #ff0000;">Find the size in Gb of the largest file</span>
3.  <span style="color: #ff0000;">Look at the first and the last 20 lines of that file.</span>

<span style="color: #ff0000;"> </span>



#### Paging through files – more and  less


Using either of these commands (the kids like less these days so use that) will page through the file and you can navigate using the space bar



There are various commands you can type while you are paging through the file

*   /asdf   will search for the next instance of the characters asdf
*   n will repeat the previous search
*   N will search backwards
*   up and down arrow will move up and down in the file onw line at a time
*   ctrl-u and ctrl-d will move up and down in the file one page at a time
*   ctrl-g will print out stats of the cursort position including line number and total lines in the file.
*   h will show the help page
*   less -S will not wrap long lines (can be useful)

There are many other options (most of which I don’t use)  Use the man less command to find out about them.



<span style="color: #ff0000;">Hands on :</span>

1.  <span style="color: #ff0000;">Using less page through the AF1.bed file in the Data directory.  Get used to the navigation commands.</span>
2.  <span style="color: #ff0000;">Find the line numbers of the first and last chr20 entries</span>





#### Odds and Ends (but very useful)

*   history          shows your command history
*   wc   <file>    returns the number of characters,words and lines in a file
*   wc -l  <file>  just returns the number of lines in a file

### 3.  That’s pretty much the basics





These commands will be most of what you need for day to day work – for example here is my recent unix history grouped by the number of times I’ve used the various commands :

10 scancel (slurm)  
10 seq  (loop stuff again)  
11 df  
11 mkdir  
11 mv  
11 wc  
12 python  
13 sort  
13 sudo   (rootytooty)  
13 tail  
15 cp  
16 find  
18 perl  
21 history  
21 sshare   (slurm)  
22 xargs  
24 rm  
26 awk  
33 for    (loops – we’ll come to these)  
40 ./iggytools/bin/seqprep  
70 vi   (my editor of choice – I’m old school)  
77 more  
85 squeue   (slurm)  
96 cd  
104 grep  
292 ls

### 4.  Getting crazy with the power.

Ok so you can log in, look around and rootle about in files and directories.      In the next session you’ll be running analysis and creating output files.   There are some extra commands that will make your life much easier when you do this and you can start to harness the power of unix.





#### 4.1 Searching using grep

Grep is a fantastic command for searching through files and directories.   The basic syntax is :

<pre>grep <pattern>  <file></pre>

So to find all entries for chr20 in our AF1.bed file we’d do

<pre>grep ‘chr20’ /n/regal/informatics/workshops/Basic_UNIX/Data/AF1.bed</pre>

Useful options   :

   -v   search for everything but the pattern

   -n   show the line number of the line found

   -C <n>  Show <n> lines context

   -r   search recursively down the directory tree (<file> is directory here)

   -i   ignore case

   -H print the filename along with the file found

   -f   find files only

   -d  find directories only

   -l   only print the filename and not the line found  (useful when there are multiple matches per file)

Extra :

Information about more complicated searching using regular expressions and egrep can be found here:

[http://ryanstutorials.net/linuxtutorial/grep.php](http://ryanstutorials.net/linuxtutorial/grep.php)

<span style="color: #ff0000;">Hands on:</span>

<span style="color: #ff0000;"> </span>

1.  <span style="color: #ff0000;">Find  how many bed peaks there are in the AF2.bed file for chr10</span>
2.  <span style="color: #ff0000;">Find which files contain the string AAAAAAAA in the Data directory</span>





#### 4.2 find



Find is a more sophisticated recursive directory search which can locate files based on a pattern but also execute commands on files when you find them.

Basic syntax is :

<pre>find <dir> -name <name pattern></pre>

This will just list all found files

Example :

<pre>find . -name “*.fastq"</pre>

The power comes when you add on the -exec option and append a command.

Example:

<pre>find . -name “*.fastq” -exec ls -l {} \;</pre>

(You put {} where the filename should go and the command should always end with \;)

This command files all fastq files and lists them

Another usage for find is finding files that are newer or older than a certain time

<pre>find . -mtime  7</pre>

Will find all files below the current directory modified more than 7 days ago

<span style="color: #ff0000;">Hands on:</span>

<span style="color: #ff0000;"> </span>

1.  <span style="color: #ff0000;">Find all files in your home directory older than 1 day</span>
2.  <span style="color: #ff0000;">Find all files under the /n/regal/informatics/workshops/Basic_UNIX directory that contain the string chr14</span>

Extra:  A good set of find examples is here [http://alvinalexander.com/unix/edu/examples/find.shtml](http://alvinalexander.com/unix/edu/examples/find.shtml)



#### 4.3 Pipes

You can feed the output of one command to the input of another using the | character

<pre><cmd1>  |  <cmd2> | <cmd3></pre>

This is best illustrated with an example

<pre>find . -name “*.bed”  |  wc -l</pre>

This will find how many *.bed files there are below the current directory

<pre>find . -name “*.bed” -exec cat {} \;  | wc -l</pre>

This finds how many lines total in all bed files found

<pre>find . -name “*.bed” -exec cat {} \; |grep chr20 | wc -l</pre>

This finds how many chr20 lines there are in all bed files



<span style="color: #ff0000;">Hands on:</span>

1.  <span style="color: #ff0000;">Find how many chr14 lines there in all bed files</span>
2.  <span style="color: #ff0000;">Find how many chr1 entries there are in all bed files under  /n/regal/informatics/workshops/Basic_Unix </span><span style="color: #ff0000;">(Hint: use the -P option and \t for a tab character)</span>



### 4.4 sort

When we have columns of data we often want to sort on a column to find the highest or lowest  entry.    A typical command looks like

<pre>sort -nk4   <file>  |less</pre>

*   n    – sort numerically
*   k4  – sort starting on the 4th column
*   k4,5  – sort using the 4th and 5th columns
*   r  – reverse sort
*   u  – sort and report unique lines

*   <span style="color: #444444; font-family: 'Open Sans';">t”,”  –  set the field delimiter to ,</span>

Example

<pre>sort -nk2 AF1.bed             # sorts the file by the 2nd column</pre>

We can have multiple column options

<pre>sort -k1,1 -k2,2n AF1.bed    # sorts the file by chromosome first and then start coord</pre>

<span style="color: #ff0000;">Hands on:</span>



1.  <span style="color: #ff0000;">Find the highest scoring 10 entries in the AF1.bed file (score is the 5th column)</span>
2.  <span style="color: #ff0000;">How many duplicate entries are there between the AF1.bed and AF2.bed files.</span>





Extra:

More sort examples are at  [http://www.theunixschool.com/2012/08/linux-sort-command-examples.html](http://www.theunixschool.com/2012/08/linux-sort-command-examples.html)





#### 4.5 awk



<span style="color: #444444; font-family: 'Open Sans', Helvetica, Arial, sans-serif;">This is where things really get powerful.  Awk is a ‘pattern scanning and processing’ utility.  It lets you search and filter files based on column and pattern.  The basic syntax is</span>

<pre>awk ‘ pattern  { action }’   filename</pre>

Again this is best shown by example

<pre>awk ‘  $1 == “chr1” { print $1,$2,$3}’   AF1.bed</pre>

Here the pattern is $1 == “chr1”  or column1 = chr1  and the action is print columns 1,2 and 3

We can do more interesting things with the pattern e.g.

<pre>awk ‘$2 > 1000000 && $3 < 2000000 { print $0}’ AF1.bed</pre>

<span style="color: #444444; font-family: 'Open Sans', Helvetica, Arial, sans-serif;">This only prints lines where the region is between 1 and 2Mb.  The $0 represents the whole line </span>

We could also do this by missing out the whole action

<pre>awk ‘$2 > 1000000 && $3 < 2000000’ AF1.bed</pre>



We can also search for substrings using the ~ character



<pre>awk ‘$1 ~ /1/' { print $1}’ AF1.bed</pre>

only prints out lines where the first column contains a 1

Similarly

<pre>awk ‘$1 !~ /1/‘ {print $1}’</pre>

<span style="color: #444444; font-family: 'Open Sans';">This prints the first column where the first column doesnt’ contain a 1</span>

There are other things you can include

<pre>awk ‘$1 ~ /^1/‘ {print $2}’</pre>

This prints the 2nd column where the 1st column starts with a 1



You can reference the line number and number of fields using NR and NF.

<pre>>awk ‘NR > 1’                      # will only print out rows 2 - end</pre>

<pre>awk ‘NF == 12’                  # will only print out rows with exactly 12 fields</pre>

<pre>awk ‘NR % 4 == 1’            # will only print lines where the line_number / 4 has a remainder of 1 (modulo)</pre>

#### Combining awk with sort and uniq

<span style="color: #444444; font-family: 'Open Sans';">The uniq command omits repeated lines.  It is often used to count repeated lines using the -c option</span>

<pre><span style="color: #444444; font-family: 'Open Sans';">awk ‘{print $1}’ Data/AF1.bed |uniq -c</span></pre>

<span style="color: #444444; font-family: 'Open Sans';">This counts the number of repeated chromosomes in the bed file (note this isn’t the total lines per chromosome as it’s only repeated lines).</span>

<span style="color: #444444; font-family: 'Open Sans';">There’s much more to awk but these commands will get you a long way.  Let’s now do something useful :</span>



<span style="color: #ff0000;">Hands on</span>

1.  <span style="color: #ff0000;">Find how many entries there are in the AF1.bed file for each chromosome</span>
2.  <span style="color: #ff0000;">Find how many entries scoring > 300 there are in the AF1.bed file for each chromosome</span>
3.  <span style="color: #ff0000;">Find the read lengths in the fastq file (Hint:  each fastq entry has 4 lines and the read length is on the 1st line of the entry)</span>

#### 4.6 perl

<span style="color: #444444; font-family: 'Open Sans', Helvetica, Arial, sans-serif;">Perl is a complete programming language but there is one aspect of it that is really useful on the command line.   This is it’s ability to search and replace strings on the fly.</span>

Syntax is

<pre>perl -pe ’s/sauage/melon/g‘   myfile</pre>

This will take each line of myfile and change each instance (the final g) of sausage into melon.   More powerfully we can have wildcards

<pre>perl -pe ’s/.*Length=//‘  myfile</pre>

This takes everything before Length= and deletes it

<span style="color: #444444; font-family: 'Open Sans', Helvetica, Arial, sans-serif;">Let’s do something more complicated</span>

<pre>perl -pe ’s/.*:(.*?)/$1/‘  myfile</pre>

This removes everything before the colon and replaces it with the string after the colon (the ? means the shortest match)

Matches can be specified more precisely than . and.*

We can use

*   \S+  for one or more  non-whitespace  (* is 0 or more and + is 1 or more)
*   \s+  for whitespace
*   \d+   for digits
*   \w+  for words
*   ^      Start of line
*   $      end of line
*   \n    newline char
*   \t    tab char
*   \.    Will match a . character (in general \ will escape a character)

This will reverse the first and 2nd words

<pre>perl -pe ’s/^(\w+)\s+(\w+)/$2 $1/‘  myfile</pre>

And of course we can use pipes

<pre>awk ‘$5 >  200 { print $1}’  |perl -pe ’s/chr//‘</pre>

This will print out all the column 1 chromosome labels for scores > 200 and strip of the chr string

<span style="color: #ff0000;">Hands on</span>

<span style="color: #ff0000;">1\.  Convert a fastq file to fasta format using one line.</span>

<span style="color: #ff0000;"> </span>

<span style="color: #ff0000;">Note: Fasta format looks like</span>

<span style="color: #ff0000;"> </span>

<span style="color: #ff0000;">>header line here</span>

<span style="color: #ff0000;">GGTTATTAGGGTGGCAGAGCCAGGAAATTGCGT</span>

<span style="color: #ff0000;">>another header line</span>

<span style="color: #ff0000;">CCTCTAAGGCGGGCCACTGTGCCAAATTCTCTA</span>

<span style="color: #ff0000;"> </span>

<span style="color: #ff0000;">2\. Extract the index strings from the Ctrl1_atac_SP2.R1.fastq fastq file and estimate the frequency distribution</span>



[http://www.softpanorama.org/Scripting/Perlorama/perl_in_command_line.shtml](http://www.softpanorama.org/Scripting/Perlorama/perl_in_command_line.shtml)



</article>