FROM centos:7

EXPOSE 80
RUN yum -y install epel-release
RUN yum -y install git python3-pip httpd

RUN echo -e "ServerName ${HOSTNAME}\n" >> /etc/httpd/conf/httpd.conf

# https://github.com/pypa/pip/issues/10219
ENV LC_ALL=en_US.UTF-8
ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt

WORKDIR /var/www
RUN git clone --single-branch --recursive https://github.com/getpelican/pelican-plugins.git && (cd pelican-plugins/jinja2content && git checkout 483215d) && mkdir website
ADD . website
RUN cd website && pelican -D content -t informatics-theme -o /var/www/html
RUN cp -r website/static/* /var/www/html

CMD ["/usr/sbin/httpd","-DFOREGROUND"]
