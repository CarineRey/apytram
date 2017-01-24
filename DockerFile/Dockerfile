FROM debian:jessie
MAINTAINER Carine Rey carine.rey@ens-lyon.fr

RUN apt-get clean && apt-get update

RUN apt-get install -qy git \
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
                        ncurses-dev


# install bowtie 2.3
RUN apt-get install -qy libtbb-dev unzip
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.0/bowtie2-2.3.0-source.zip &&\
    unzip bowtie2-2.3.0-source.zip && cd bowtie2-2.3.0/ && make && make install


# install ncbi-blast+=2.6.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz &&\
    tar zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-2.6.0+/bin/:$PATH


### install openjdk-8-jdk
RUN echo "deb http://http.debian.net/debian jessie-backports main" >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -qy openjdk-8-jdk
#update-java-alternatives --set java-1.8.0-openjdk-amd64


### install Trinity
WORKDIR /opt
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz; \
     tar zxvf Trinity-v2.3.2.tar.gz; rm  Trinity-v2.3.2.tar.gz
WORKDIR /opt/trinityrnaseq-Trinity-v2.3.2/
RUN make
ENV PATH $PATH:/opt/trinityrnaseq-Trinity-v2.3.2/:/opt/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/

#### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram
ENV PATH $PATH:/opt/apytram


RUN apt-get install -qy vim htop

RUN mkdir /data
WORKDIR /data
