############################################################
# Dockerfile to build CRISPResso2
############################################################

#FROM continuumio/miniconda3
FROM mambaorg/micromamba:1.5.0

USER root
# File Author / Maintainer
MAINTAINER Kendell Clement
RUN apt-get update && apt-get install gcc g++ bowtie2 samtools libsys-hostname-long-perl \
  -y --no-install-recommends \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && micromamba install -c defaults -c conda-forge -c bioconda -y -n base --debug fastp "numpy<2" cython jinja2 tbb=2020.2 pyparsing=2.3.1 scipy matplotlib-base pandas plotly\
  && micromamba clean --all --yes

#install ms fonts
RUN echo "deb http://httpredir.debian.org/debian buster main contrib" > /etc/apt/sources.list \
  && echo "deb http://security.debian.org/ buster/updates main contrib" >> /etc/apt/sources.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# install crispresso
COPY . /CRISPResso2
WORKDIR /CRISPResso2

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python setup.py install \
  && CRISPResso -h \
  && CRISPRessoBatch -h \
  && CRISPRessoPooled -h \
  && CRISPRessoWGS -h \
  && CRISPRessoCompare -h

USER $MAMBA_USER

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python","/CRISPResso2/CRISPResso2_router.py"]