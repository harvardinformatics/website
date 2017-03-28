Title: Unix Command Line Tips and Tricks
Date: 2015-01-01 10:00
Category: Tutorials
Tags: Unix
Summary: This basic unix tutorial will get you up and running for making files, listing directories and running single commands.   It will show you what you need for about 80% of your command line work.

This basic unix tutorial will get you up and running for making files, listing directories and running single commands.   It will show you what you need for about 80% of your command line work.

If you find yourself at the command line a lot and have a lot of files to process then there are built in features on the bash command line (or shell) to help you.

The following is a list of tips and tricks that have proved useful to us and hopefully will be to you.

### 1 Files and filenames

#### 1.1  Changing the extension of a filename (and basic variable string manipulation)

To do this we want to take a filename, take off its extension and then put a fresh one on.   The syntax to do this is easy but not particularly memorable

First put the filename into a variable.   (In practice you wouldn’t do this for a single file – it would be much easier just to enter the “mv myfile.txt myfile.dat” command.  This syntax comes into play, however when we’re processing a lot of files or doing some processing using a script.   This will come later so bear with us for now.)

So – back to setting the variable i to a filename.

    :::bash
    i=myfile.txt

Check we have the right text in the variable.

    :::bash
    echo $i
    > myfile.txt

Note we set the variable with no dollar sign and reference it with a dollar sign.  I have no idea why.

Next we use a nifty bash trick to take off the file extension and put it into a variable `$j`

    :::bash
    j=${i%.*}
    echo $j
    > myfile

And finally pop the new extension on

    :::bash
    k=$j.dat
    echo $k
    
    > myfile.dat

I know this seems long winded but at long last we rename the file

    :::bash
    mv $i $k

A much shorter (but not quite so understandable) version is to combine everything together and

    :::bash
    mv $i ${i%.*}.dat

#### 1.2 Manipulating filenames (strings) on the command line

There are several other similar tricks to the one in section 1.1 to manipulate filenames (or strictly speaking strings).

    :::bash
    ${variable%pattern}  Trim the shortest match from the end
    ${variable##pattern} Trim the longest match from the beginning
    ${variable%%pattern} Trim the longest match from the end
    ${variable#pattern}  Trim the shortest match from the beginning


We’ll use some of these in the next few sections.

#### 1.3 Getting the extension of a filename

This uses a similar command to the way we changed the extension of a filename but now we use the command that matches from the beginning and not from the end

    :::bash
    i=${i##*.}
    echo $i
    > txt

#### 1.4  Getting the directory from a full path

    :::ash
    dirname=`dirname “$file”`

or

    :::bash
    dirname=${file%/*}

i.e. everything before the last ‘/’

#### 1.5  Removing the directory from a full path

Finally for completeness it is highly likely you’ll want to process filenames that have their full path attached (`/n/home_rc/mclamp/sausage.dat` for instance).  In these cases there are a couple of ways to do just extract the filename without the directory :

    :::bash
    i=${i##*/}

or using the built in basename function

    :::bash
    i=/n/home_rc/mclamp/myfile.txt
    i=`basename $i`
    echo $i
    > myfile.txt

Note here we’re using the ubiquitous back ticks \`..\` which executes whatever is inside them and replaces the command with the output.

#### 1.6  Backing up a file using brace expansion

A well loved use of brace expansion is to make a backup of a file

    :::bash
    cp myfile.txt{,.bak}

(Note: no spaces!)

#### 1.7 Making a series of files

The shell has a number of ways of generating a series of things.   First let’s start with numbers.  If we want a range of numbers 1 – 10 we enclose them in curly brackets and separate them with .. (this is usually referred to as brace expansion).

    :::bash
    echo {1..10}
    1 2 3 4 5 6 7 8 9 10

We can use this in commands  – for instance to create a set of files :

    :::bash
    touch pog{1..10}
    ls pog{1..10}
    pog1   pog10  pog2   pog3   pog4   pog5   pog6   pog7   pog8   pog9

And to tidy up…

    :::bash
    rm pog{1..10}

We can also use this to good effect in loops (of which more later)

    :::bash
    for i in 1..5 ; do
    echo $i
    done
    1
    2
    3
    4
    5


#### 1.8  Processing a set of files (for loops by stealth)

If we’re doing lots of processing then this often means iterating over a number of files – maybe a number of .fastq files output from a sequencing run perhaps. This is where we bring in loops – let’s just dive in.

    :::bash
    for i in *.fastq ; do
      j=${i%.fastq}.out
      run_something $i > $j
    done


This takes all the .fastq files in the directory, changes the extension to .out, runs a command and puts the output in the new file

#### 1.9 Processing a more complicated list of files

Instead of just a straight vanilla wildcard filename expansion *.fastq we can get quite exotic and put commands, pipes, greps in the list part.

Example :

    :::bash
    for i in `find . -name “*.fastq” |grep R1 ` ; do    # Note the backticks!
    
       do_something_with_R1_fastq_files_here
    
    done


This looks for all .fastq files but filters them using grep for only those containing the characters R1

### 2 Manipulating the contents of files

#### 2.1 Even more complicated loops

We’re not limited to filenames in the for loop list. We can do fancy things with the contents of files

    :::bash
    for i in `cat *.dat | awk ‘$1 == “SEQ” { print $5}’` ; do
    
        echo $i
    
    done


So here we’re looking through the contents of all files ending with .dat (the cat *.dat) and filtering it through awk. If the first field in the file is SEQ then print out the 5th column.

This is actually more useful than you might think. Let’s do something with fastq files again

    :::bash
    for i in `cat *.fastq|awk ‘NR%4 == 1’ ` ; do
    
       # Hmm maybe not a good example here
    
    done


#### 2.2 Once more with feeling – using awk

We leapt ahead a bit there. Loops, backticks, awk, grep etc all in one thing. Let’s look at awk a bit more first.

awk is at heart a way of filtering text. Many programs produce text output in columns and awk is an excellent way to filter and process the output

A very basic (yet still useful) use of awk is to output columns from a file. e.g.

    :::bash
    awk '{ print $3, $4, $10 }' myfile.dat

This will print columns 3,4 and 10 only from myfile.dat (remember the curly brackets folks!)

If we put something before the curly brackets this acts as an ‘if’ statement. For instance

    :::bash
    awk ' $2 == "SEQ" { print $3,$4,$10}' myfile.dat

This only prints out columns 3,4 and 10 if SEQ is in column 2\. Similarly we can put numerical comparisons here too

    :::bash
    awk ' $1 < 0.05 { print $6}' myfile.dat

This filters the file to rows where column 1 < .05 and then prints column 6

If we want to get fancier let's combine the filename manipulation with awk to pipe into another file

    :::bash
    i=myfile.dat
    awk ' $1 < 0.05 { print $0 }' $i > ${i%.*}.under05.dat


Here $0 will print the whole line and we end up with the output in myfile.under05.dat. Nice and neat yes?

We can also search for partial string matches using ~ /mystring/

    :::bash
    awk ' $2 ~ /metal_ion_binding/  {print $0 }' $i > ${i%.*}.metal_ion_binding.dat

#### 2.3 Fun with awk - summing and averaging a column

As well as the 'if' section and the 'main' section of an awk statement we can do things before and after filtering. For instance if we want the sum of a column in a file

    :::bash
    awk '{ s += $4} END {print s}'
    471.948


Here s is being used as an awk variable and keeps the total of column 4\. At the end of the file we print out the total

Similarly we can do stuff at the beginning

    :::bash
    awk 'BEGIN {print "Total"} { s += $4} END {print s}'
    Total
    471.948


And just to be a little fancy we can combine all sorts of things

    :::bash
    awk 'NR % 4 == 2 { s += length($1); t++} END {print s/t}'

Explanation NR - row number, length($1) returns the length of the string.  
So on rows 2,6,10 etc we sum the length of column 1 and keep a tally of entries summed in variable t. At the end we print the average of the length of column 1

This is an actual command to calculate the average length of a sequencing read in a fastq file

    :::bash
    awk 'NR % 4 == 2 { s += length($1); t++} END {print s/t}'  ../pogpipe/testdata/sample_1.fq
    36


(Admittedly this is a pretty old file these days - 36 base reads? Pffft)

### 3. Date strings

#### 3.1 Using datestamps in filenames

If I’m running lots of things I quite often want to put output in datestamped files (for example putting output into a new file every day). Using command substitution we can do this quite easily.

The date command takes a string with codes for day, month, year etc.

    :::bash
    today=`date +%d-%b-%Y`
    echo $today
    
    13-Jun-2014


We can also put this directly into a filename :

    :::bash
    ls -ltra /n/home_rc/mclamp/ > /tmp/dirlist.$(date +%d-%b-%Y).log</pre>

which creates

    :::bash
    /tmp/dirlist.13-Jun-2014.log</pre>

Note here we could have used back ticks but we used the $(command) syntax instead. Either can be used but I prefer `..` which probably means the favored way is $(..) [[Edit: the favored way is indeed $(command) and for the reason that the opening and closing characters are different so it makes the code easier to read. I have to grudgingly admit they have a point.]]

#### 3.2 Other date formats

If we want a purely numeric date string YYYY-mm-dd use

    :::bash
    today=$(date +%Y-%m-%d)
    echo $today
    
    2014-06-13


And if we want hours:mins:secs we do

    :::bash
    now=$(date +Y-%m-%d_%H:%M:%S)
    echo $now
    
    2014-06-13_16:36:23


A couple of extras with shortcuts for common formats

    :::bash
    date +%Y-%m-%d\ %R

gives

    :::bash
    2011-03-25 09:48

    :::bash
    date +%Y-%m-%d\ %T 

gives

    :::bash
    2011-03-25 09:51:05

Now you’ll never have an unstamped log file ever again.

### 4. Miscellaneous useful stuff

#### 4.1 Fun with loops - turn a command into a continuous monitor

Take any command and turn it into a monitor in a terminal window

Keep an eye on who is logging in

    :::bash
    while [ 1 ]; do clear; last |head -10; sleep 3; done

Keep an eye on who is logged in

    :::bash
    while [ 1 ]; do clear; w; sleep 3; done

#### 4.2 Bash command line navigation

This is mostly lifted from http://splike.com/wiki/Bash_Scripting_FAQ

M is the ‘meta’ key - on my macbook it is the escape key but will be the alt key on linux keyboards (is this true - I don’t have one to hand). Also on a macbook you have to release the meta key before pressing the next one. It’s awkward to start but gets easier over time.

    :::bash
    M-f     forward one word
    M-b     back one word
    M-d     delete one word

    so to delete the word you're on is M-b M-d i.e. keep you finger on the alt key and then press b followed by d.

    or!!

    M-C-h  backward kill word

    M-<   beginning of history

    C-r reverse search history
    C-s forward search history

    C-M-y yank first argument from previous line.   The man says it can take n as an argument but I can't make this work

    M-.  or M-_   yank last arg

    C-t  transpose chars (previous char and current char)
    M-t  transpose words

    M-u upcase word
    M-l lowercase word

    M-c capitalize word

    C-k delete from cursor to end of line (k for kill)

    M-d kill word
    M-del backward kill word
    C-x-del  delete from cursor to beginning of line


#### 4.3 Reading input from the command line

    :::bash
    read -p “Enter your name (first last) : “ first name last name
    
    Enter your name (first last) : Michele Clamp


You can also read input from stdin in a loop

    :::bash
    for i in read tmp  ; do
    
      echo $tmp
    
    done

