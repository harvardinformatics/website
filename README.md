# website
Website markdown files


# Install pelican and markdown (python 2.7 required)

pip install pelican markdown

# Install the pelican plugins (the pelicanconf.py file assumes it is cloned at the same level as this repo)

cd ..
git clone --recursive https://github.com/getpelican/pelican-plugins

#From the root of the project run the conversion with the specified theme:

cd website 
pelican content -t `pwd`/informatics-theme

#Then go to the output directory and start the server

cd output && python -m pelican.server

#Should be visible from localhost:8000

