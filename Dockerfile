FROM centos:7 AS dev

EXPOSE 80
RUN yum -y install epel-release
RUN yum -y install git python python2-pip httpd

RUN echo -e "ServerName ${HOSTNAME}\n" >> /etc/httpd/conf/httpd.conf

ADD requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt

WORKDIR /var/www
# work around python2 compatibility issue in jinja2content.py
# support tipue.js 5.x bundled in informatics-theme
RUN git clone https://github.com/getpelican/pelican-plugins.git \
 && cd pelican-plugins \
 && git checkout 736814b \
 && git show 130bdd8:jinja2content/jinja2content.py > jinja2content/jinja2content.py \
 && git show 2dcdca8:tipue_search/tipue_search.py > tipue_search/tipue_search.py \
 && git submodule update --init -- pelican-toc rmd_reader pin_to_top \
 && rm -rf .git

FROM dev AS final

ADD . website
RUN cd website && pelican -D content -t informatics-theme -o /var/www/html
RUN cp -r website/static/* /var/www/html

CMD ["/usr/sbin/httpd","-DFOREGROUND"]
