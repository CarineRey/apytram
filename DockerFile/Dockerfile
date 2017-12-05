FROM debian:stretch-slim
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

### install openjdk-8-jre
RUN echo "deb http://httpredir.debian.org/debian stretch-backports main" >> /etc/apt/sources.list && mkdir /usr/share/man/man1 && \
    apt-get update && apt-get install --no-install-recommends -t stretch-backports -qy ca-certificates-java openjdk-8-jre-headless 


###  install bowtie2 2.3.3.1
ENV bowtie_version 2.3.3.1
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/"$bowtie_version"/bowtie2-"$bowtie_version"-source.zip &&\
    unzip bowtie2-"$bowtie_version"-source.zip && cd bowtie2-"$bowtie_version"/ && make && make install


###  install ncbi-blast+=2.6.0
ENV blast_version 2.6.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/"$blast_version"/ncbi-blast-"$blast_version"+-x64-linux.tar.gz &&\
    tar zxf ncbi-blast-"$blast_version"+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-"$blast_version"+/bin/:$PATH


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
ENV transdecoder_version 3.0.1
WORKDIR /opt
RUN wget https://github.com/TransDecoder/TransDecoder/archive/v"$transdecoder_version".zip &&\
    unzip v"$transdecoder_version".zip  && rm v"$transdecoder_version".zip && \
    ls /opt/ && cd /opt/TransDecoder-"$transdecoder_version"/ && make && ls /opt
ENV PATH /opt/TransDecoder-"$transdecoder_version"/:$PATH

## install R and ggplot2 for TransDecoder
#RUN echo 'install.packages(c("ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE); \
#    source("https://bioconductor.org/biocLite.R");\
#    biocLite("seqLogo")' > /tmp/install.R && Rscript /tmp/install.R && rm /tmp/install.R

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

