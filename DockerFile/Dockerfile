FROM debian:jessie
MAINTAINER Carine Rey carine.rey@ens-lyon.fr

RUN apt-get clean && apt-get update

RUN apt-get install --no-install-recommends  -qy git \
                        wget \
                        build-essential \
                        mafft \
                        exonerate \
                        python-numpy \
                        python-pandas \
                        python-matplotlib \
# to install trinity
                        samtools=0.1.19-1 \
                        zlib1g-dev \
                        ncurses-dev \
# to install bowtie2
                        libtbb-dev \
                        unzip \
# to user
                        vim \
                        htop


### install openjdk-8-jre
RUN echo "deb http://httpredir.debian.org/debian jessie-backports main" >> /etc/apt/sources.list
RUN apt-get update && apt-get install --no-install-recommends -t jessie-backports -qy ca-certificates-java openjdk-8-jre-headless 

    
###  install bowtie 2.3
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.0/bowtie2-2.3.0-source.zip &&\
    unzip bowtie2-2.3.0-source.zip && cd bowtie2-2.3.0/ && make && make install


###  install ncbi-blast+=2.6.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz &&\
    tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-2.6.0+/bin/:$PATH


### install Trinity
WORKDIR /opt
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz && \
     tar zxf Trinity-v2.3.2.tar.gz && rm Trinity-v2.3.2.tar.gz
WORKDIR /opt/trinityrnaseq-Trinity-v2.3.2/
RUN make
ENV PATH $PATH:/opt/trinityrnaseq-Trinity-v2.3.2/:/opt/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/

### install NextGenMap
RUN git clone https://github.com/Cibiv/NextGenMap.git  /opt/NextGenMap/
RUN mkdir /opt/NextGenMap/build
WORKDIR /opt/NextGenMap/build
RUN apt-get install --no-install-recommends  -qy cmake
RUN cmake .. && make
RUN make install prefix=/usr/local/bin
ENV PATH /opt/NextGenMap/bin/ngm-0.5.3:$PATH


### install seqtk
RUN git clone https://github.com/lh3/seqtk.git  /opt/seqtk/
WORKDIR /opt/seqtk
RUN make
ENV PATH /opt/seqtk:$PATH


### install cd-hit
RUN git clone https://github.com/weizhongli/cdhit.git /opt/cdhit/
WORKDIR /opt/cdhit/cd-hit-auxtools
RUN make
ENV PATH /opt/cdhit/cd-hit-auxtools:$PATH
WORKDIR /opt/cdhit
RUN make
RUN make install


### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram
ENV PATH /opt/apytram:$PATH


### 
RUN mkdir /data
WORKDIR /data
