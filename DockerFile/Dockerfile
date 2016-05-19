FROM debian:jessie
MAINTAINER Carine Rey carine.rey@ens-lyon.fr
RUN apt-get update && \
    apt-get install -qy git \
                        wget \
                        build-essential \
						mafft \
						exonerate \
						python-numpy \				
					    python-pandas \
	                    python-matplotlib \
# to install trinity	                    
	                    bowtie \ 
	                    openjdk-7-jdk \
	                    samtools=0.1.19-1 \ 
	                    zlib1g-dev \
	                    ncurses-dev 

# install ncbi-blast+=2.2.28-2
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz &&\
    tar zxvf ncbi-blast-2.2.28+-x64-linux.tar.gz 
ENV PATH /opt/ncbi-blast-2.2.28+/bin/:$PATH

### install Trinity
WORKDIR /opt
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz; \
     tar zxvf v2.1.1.tar.gz; rm  v2.1.1.tar.gz 
WORKDIR /opt/trinityrnaseq-2.1.1/
RUN make

#### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram     
ENV PATH /opt/apytram:/opt/trinityrnaseq-2.1.1/:/opt/trinityrnaseq-2.1.1/trinity-plugins/:$PATH


RUN apt-get install -qy vim htop

RUN mkdir /data
WORKDIR /data
