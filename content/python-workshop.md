Title: Odyssey Python Workshop
Date: 2016-1-13 1:00
Category: Tutorials
Author: aaron_kitzmiller@harvard.edu
Tags: Linux, Workshop, Python
Summary: This tutorial is intended for people who are familiar with the basics of unix and Odyssey and want to know how to use Python on the cluster.




## 1. Introduction and prerequisites

This tutorial is intended for people who are familiar with the basics of unix and Odyssey and want to know how to use Python on the cluster.

*   A Harvard FAS RC cluster account 
    [https://portal.rc.fas.harvard.edu/request/account/new/](account-request>)

*   A terminal program so you can log into the cluster 
    [https://rc.fas.harvard.edu/resources/access-and-login/](access-and-login>)

*   The ability to log into the cluster 
    [https://rc.fas.harvard.edu/resources/access-and-login/](access-and-login>)

*   <b>Familiarity with basic unix commands involving files (cp, rm , mv, cat, less, ls)
and/or have attended the ‘Basic Unix’ workshop previously.</b>

*   Basic skills with a Linux editor such as vi, nano, emacs, etc.



## 2. Summary of Topics covered

## 3. Python background
Short background about interpreted languages, philosophy of Python, impact on scientific computing, packages and pypi, python2 versus python3, REPL


    :::python
    >> # Hello world with a loop
    >> for i in xrange(1,10):
    >>    print "Hello, World, my name is %s." % 'Aaron Kitzmiller'
    >> print "Done."

    :::python
    >> # Hello world with a formatted string
    >> for i in xrange(1,10):
    >>    print "Hello, World, my name is %s for the %dth time." % ('Aaron Kitzmiller', i)
    >> print "Done."

    (Show the string formatting page here, so they know where to look.)

    :::python
    >> # Hellow world with bad indent
    >> for i in xrange(1,10):
    >>     print "Hello"
    >>   print "World."
    >> print "Done"

 



