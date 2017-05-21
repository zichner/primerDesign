# use the ubuntu base image
FROM ubuntu:16.04

MAINTAINER Thomas Zichner


# install required packages
RUN apt-get update && apt-get install --no-install-recommends -y \
	primer3 \
	ncbi-blast+ \
    && apt-get clean \
	&& apt-get purge

RUN apt-get install --no-install-recommends -y python-pip python-setuptools \
	&& pip install --upgrade pip \
	&& pip install xlrd xlwt \
	&& apt-get autoremove -y python-pip python-setuptools

RUN apt-get install --no-install-recommends -y aria2 \
	&& mkdir /opt/blastDb \
	&& cd /opt/blastDb \
	&& aria2c -x 8 -s 8 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz \
	&& gunzip hs37d5.fa.gz \
	&& makeblastdb -in hs37d5.fa -dbtype nucl -parse_seqids \
	&& rm hs37d5.fa \
	&& apt-get autoremove -y aria2

RUN apt-get install --no-install-recommends -y git \
	&& cd /opt \
	&& git clone https://github.com/zichner/primerDesign.git \
	&& apt-get autoremove -y git

RUN mkdir /data /blastDb

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# Change workdir to /data/
WORKDIR /data/

# by default, primerDesign is executed
ENTRYPOINT ["python", "/opt/primerDesign/primerDesign.py"]

