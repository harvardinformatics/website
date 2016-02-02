Title: Port Forwarding on Odyssey Nodes
Date: 2015-01-01
Author: Aaron Kitzmiller
Category: Tutorials
Tags: screencast, Odyssey, IPython
Summary:  Using the Odyssey slurm tunneling plugin to simplify port-forwarding

Harvard Informatics has developed a plugin for the Odyssey Slurm system that simplifies port forwarding between submit and execution hosts.  As described below, this can be particularly useful for iPython notebooks.

Communication between computers is typically done through ports, numbered memory locations that listen for TCP (and UDP) traffic. Some of the most commonly used ports are 22 (for ssh) and 80 (http). Port forwarding, or tunneling, is a technique that allows communication on different ports (e.g. 8888) to be passed over an ssh connection (port 22). If, for example, ssh connections were allowed between two computers, but not web traffic, you could still run a web server and tunnel the port 80 traffic over an ssh connection. 

<figure>
	    <a class="img" href="/images/tunnel-example.jpg">
		    <img class="img-responsive" src="/images/tunnel-example.jpg"></img>
	    </a>
<figcaption>Example of port forwarding of web traffic over an ssh connection when web traffic is blocked</figcaption>
</figure>   

A new srun parameter, `--tunnel`, allows a user to specify one or more pairs of port numbers that will be tunneled through a standard ssh connection. 

    :::bash
    srun --pty -p interact --mem 4000 -t 300 --tunnel 8889:8888 bash
    
This will start a new interactive session in which port 8889 on the submit host will be listening to port 8888 on the execution host.  When you exit the interactive session, the tunnel will be terminated. 

If you need more than one port forwarded, you can provide a comma separated list 

    :::bash
    srun --pty -p interact --mem 4000 -t 300 --tunnel 8889:8888,8000:8000 bash 
    
There are some limitations to the plugin.  You can only launch a single srun command that includes a --tunnel parameter on a given submit host.  Again, since you can provide a list of port pairs, this shouldn't be necessary.

## Setting up an iPython notebook using port forwarding

Running an iPython notebook is a common scientific use case for port forwarding. iPython notebooks are typically viewed by pointing a browser at "localhost:8888" after running the command `ipython notebook`. 

On a shared system like Odyssey, notebooks should be executed on compute hosts rather than login nodes. While it is possible with x11 forwarding to run a notebook and browser on the compute host, the performance of x11 forwarding can be poor. It is much more efficient to forward the web server traffic from the compute host back to the login node and then to your desktop, or from the compute host back to an NX desktop. 

This screencast demonstrates how to use the port forwarding plugin from an [NX desktop on Odyssey](https://rc.fas.harvard.edu/resources/access-and-login/#Consider_an_NX_remote_desktop_for_graphical_applications_like_Matlab_and_RStudio "NX desktops on Odyssey"). 

[embed]http://youtu.be/dVXJRm6GLY4[/embed] 

If instead, you'd like to connect directly to your laptop / desktop, this screencast will show you how. 

[embed]http://youtu.be/ITy-Ow77A5g?list=UUDslFRDz5P9-069ktbuEkRQ[/embed]