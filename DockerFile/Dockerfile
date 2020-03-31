FROM debian:buster-slim
MAINTAINER Carine Rey carine.rey@ens-lyon.org

RUN apt-get update

RUN apt-get install --no-install-recommends  -qy git \
                        wget \
                        cmake \
                        build-essential \
                        mafft \
                        exonerate \
                        python-pip \
                        python-dev \
                        python-setuptools \
#                       python-numpy \
#                       python-pandas \
#                       python-matplotlib \
#                       python-biopython \
# to install trinity
                        zlib1g-dev \
                        ncurses-dev \
                        rsync \
                        salmon \
                        jellyfish \
                        samtools\
#To install Transdecoder
                        libany-uri-escape-perl \
#To install R
#                        r-base \
#                        curl \
#                        libcairo2-dev \
#                        libssl-dev \
#                        libcurl4-openssl-dev \
# to install bowtie2
                        libtbb-dev \
                        unzip \
# to user
                        vim \
                        htop \
                        xvfb \
                        gosu 


RUN pip install numpy
RUN pip install pandas
RUN pip install matplotlib
RUN pip install biopython

### install openjdk-11-jre
RUN echo "deb http://httpredir.debian.org/debian stretch-backports main" >> /etc/apt/sources.list && mkdir /usr/share/man/man1 && \
    apt-get update && apt-get install --no-install-recommends -t stretch-backports -qy ca-certificates-java openjdk-11-jre-headless 


###  install bowtie2 2.3.5.1
ENV bowtie_version 2.3.5.1
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/"$bowtie_version"/bowtie2-"$bowtie_version"-source.zip &&\
    unzip bowtie2-"$bowtie_version"-source.zip && cd bowtie2-"$bowtie_version"/ && make -j 5 && make install


###  install ncbi-blast+=2.9.0
ENV blast_version 2.9.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/"$blast_version"/ncbi-blast-"$blast_version"+-x64-linux.tar.gz &&\
    tar zxf ncbi-blast-"$blast_version"+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-"$blast_version"+/bin/:$PATH


### install Trinity
ENV trinity_version v2.9.1
WORKDIR /opt

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/"$trinity_version"/trinityrnaseq-"$trinity_version".FULL.tar.gz && \
     tar zxf trinityrnaseq-"$trinity_version".FULL.tar.gz && rm trinityrnaseq-"$trinity_version".FULL.tar.gz &&\
     cd /opt/trinityrnaseq-"$trinity_version"/ && make -j 15 && make install
ENV PATH $PATH:/opt/trinityrnaseq-"$trinity_version"/:/opt/trinityrnaseq-"$trinity_version"/trinity-plugins/


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
ENV transdecoder_version TransDecoder-v5.5.0
WORKDIR /opt

RUN wget https://github.com/TransDecoder/TransDecoder/archive/"$transdecoder_version".zip &&\
    unzip "$transdecoder_version".zip  && rm "$transdecoder_version".zip && \
    ls /opt/ && cd /opt/TransDecoder-"$transdecoder_version"/ && make && ls /opt
ENV PATH /opt/TransDecoder-"$transdecoder_version"/:$PATH


### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram
ENV PATH /opt/apytram:$PATH
RUN cd /opt/apytram && make test

##### install entrypoint
ADD entrypoint.sh /usr/local/bin/entrypoint.sh

##### build a working directory
WORKDIR /data

CMD ["apytram.py", "-h"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

