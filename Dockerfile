FROM centos:7

RUN yum update && yum -y install git python python-pip httpd
ADD requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt


