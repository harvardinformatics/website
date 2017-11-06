# Harvard FAS Informatics website
Website markdown files.  

## Editing the website is easy
Just checkout the repository, edit a markdown file and check it back in.  Which, assuming you know Markdown and all of the conventions for the website, is pretty easy.  

If you're going to make some significant changes, beyond adding a page or modifying an existing one, you should create a branch and then do a merge request so that someone else can take a look.

## Building the website the easy way
So, now, the easiest way to build the website is to use Docker.  

### Install Docker
This is pretty easy if you have a Mac.  The official DMG is pretty good.  You shouldn't need to homebrew.

### Checkout the repository and enter the top level directory
    git clone git@github.com:harvardinformatics/website.git
    cd website

### Build a docker image
The "tag" that you give it can be anything, but "website" makes sense.  The Docker build will take care of the setup of the Python environment, getting the pelican-plugins, and converting the source Markdown files into HTML.

    docker build -t website .
    
### Run it on your local machine
The container includes a web server listening to port 80.  If you run the container and map a port on your localhost, you can view it there. To see it at http://localhost:8080, 

    docker run -d -p 8080:80 website

Regardless of whether you are using Docker, once you're satisfied, just check your changes in (or merge your branch into master).  The live website should be automatically updated.

## Building the website the hard way
If you want to build it without Docker, you have to do all this setup stuff.

### Install pelican and markdown (python 2.7 required)

    pip install pelican markdown beautifulsoup4 icalendar

### Install the pelican plugins (the pelicanconf.py file assumes it is cloned at the same level as this repo)

    cd ..
    git clone --recursive https://github.com/getpelican/pelican-plugins

If clone --recursive does not work (you might see a complaint about a non-existent plugin), you may need to manually init and update the submodule:

    cd pelican-plugins/pelican-toc
    git submodule init
    git submodule update

### From the root of the project run the conversion with the specified theme:

    cd website 
    pelican content -t `pwd`/informatics-theme -o output

### Then go to the output directory and start the server

    cd output && python -m pelican.server

### Or use a one-liner
    pelican content -t `pwd`/informatics-theme -o output && (cd output && python -m pelican.server)

### Should be visible from localhost:8000


### Making live on the website

rsync is no longer necessary to make your changes live.  If you push to the master branch, the website is automatically updated by a github hook.  You pretty much have no excuse.

