FROM debian:jessie-slim
MAINTAINER Carine Rey carine.rey@ens-lyon.org

RUN apt-get update

RUN apt-get install --no-install-recommends  -qy git \
                        wget \
                        cmake \
                        build-essential \
                        mafft \
                        exonerate \
                        python-numpy \
                        python-pandas \
                        python-matplotlib \
                        python-biopython \
# to install trinity
                        zlib1g-dev \
                        ncurses-dev \
                        rsync \
#To install Transdecoder
                        libany-uri-escape-perl \
# to install bowtie2
                        libtbb-dev \
                        unzip \
# to user
                        vim \
                        htop


### install openjdk-8-jre
RUN echo "deb http://httpredir.debian.org/debian jessie-backports main" >> /etc/apt/sources.list && mkdir /usr/share/man/man1 && \
    apt-get update && apt-get install --no-install-recommends -t jessie-backports -qy ca-certificates-java openjdk-8-jre-headless 


###  install bowtie 2.3.2
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-source.zip &&\
    unzip bowtie2-2.3.2-source.zip && cd bowtie2-2.3.2/ && make && make install


###  install ncbi-blast+=2.6.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz &&\
    tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-2.6.0+/bin/:$PATH


### install Trinity
ENV trinity_version v2.5.1
WORKDIR /opt
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-"$trinity_version".tar.gz && \
     tar zxf Trinity-"$trinity_version".tar.gz && rm Trinity-"$trinity_version".tar.gz &&\
     cd /opt/trinityrnaseq-Trinity-"$trinity_version"/ && make && make install
ENV PATH $PATH:/opt/trinityrnaseq-Trinity-"$trinity_version"/:/opt/trinityrnaseq-Trinity-"$trinity_version"/trinity-plugins/


### install seqtk
WORKDIR /opt/seqtk
RUN git clone https://github.com/lh3/seqtk.git  /opt/seqtk/ && make
ENV PATH /opt/seqtk:$PATH


### install cd-hit
RUN git clone https://github.com/weizhongli/cdhit.git /opt/cdhit/
RUN cd /opt/cdhit/cd-hit-auxtools &&  make && \
    cd /opt/cdhit &&  make  && make install
ENV PATH /opt/cdhit/cd-hit-auxtools:$PATH


#### install TransDecoder
ENV transdecoder_version TransDecoder-v5.0.2
WORKDIR /opt
RUN wget https://github.com/TransDecoder/TransDecoder/archive/"$transdecoder_version".zip &&\
    unzip "$transdecoder_version".zip  && rm "$transdecoder_version".zip && \
    ls /opt/ && cd /opt/TransDecoder-"$transdecoder_version"/ && make && ls /opt
ENV PATH /opt/TransDecoder-"$transdecoder_version"/:$PATH

### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram
ENV PATH /opt/apytram:$PATH

### 
WORKDIR /data
