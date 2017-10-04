FROM centos:7

EXPOSE 80
RUN yum -y install epel-release
RUN yum -y install git python python2-pip httpd

RUN echo -e "ServerName ${HOSTNAME}\n" >> /etc/httpd/conf/httpd.conf

ADD requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt

WORKDIR /var/www
RUN git clone --recursive https://github.com/getpelican/pelican-plugins.git && git clone https://github.com/harvardinformatics/website.git
RUN cd website && pelican content -t informatics-theme -o /var/www/html

CMD ["/usr/sbin/httpd","-DFOREGROUND"]
