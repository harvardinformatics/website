## Harvard FAS Informatics website
Website markdown files


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

### Should be visible from localhost:8000


### Making live on the website


For updating the web site, check out the files, make your edits, and do a make rsync_upload
You'll want to set SSH_USER to your username
