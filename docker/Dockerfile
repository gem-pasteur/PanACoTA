from ubuntu:18.04

# Update apt-get packages
RUN apt-get update &&\
    apt-get -y upgrade


# Install packages needed and update pip
RUN apt-get install -y \
        wget \
        python3 \
        python3-pip
# Upgrade pip
RUN pip3 install --upgrade pip


# Update makeblastdb and blastp for prokka
WORKDIR /tmp
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz &&\
    tar zxvpf ncbi-blast-2.10.1+-x64-linux.tar.gz &&\
    cp /tmp/ncbi-blast-2.10.1+/bin/makeblastdb /usr/local/bin/ &&\
    cp /tmp/ncbi-blast-2.10.1+/bin/blastp /usr/local/bin/ &&\
    rm ncbi-blast-2.10.1+-x64-linux.tar.gz


# Install mash
WORKDIR /tmp
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar &&\
    tar -xf mash-Linux64-v2.2.tar &&\
    rm mash-Linux64-v2.2.tar &&\
    mv /tmp/mash-Linux64-v2.2/mash /usr/local/bin &&\
    rm -r mash-Linux64-v2.2


# Install barrnap
WORKDIR /tmp
RUN wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz &&\
    tar -xf 0.8.tar.gz &&\
    rm 0.8.tar.gz &&\
    mv /tmp/barrnap-0.8/bin/barrnap /usr/local/bin  &&\
    # Remove heavy useless files
    rm -r /tmp/barrnap-0.8/examples /tmp/barrnap-0.8/build/*.aln

# Install prodigal
WORKDIR /tmp
RUN wget https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz &&\
    tar -xzf v2.6.3.tar.gz &&\
    rm v2.6.3.tar.gz
WORKDIR /tmp/Prodigal-2.6.3
RUN make &&\
    make install

# Install prokka:
WORKDIR /tmp
RUN apt-get install -y\
        libdatetime-perl \
        libxml-simple-perl \
        libdigest-md5-perl \
        git \
        default-jre \
        bioperl
RUN git clone https://github.com/tseemann/prokka.git
RUN /tmp/prokka/bin/prokka --setupdb &&\
    ln -s /tmp/prokka/bin/prokka /usr/local/bin


# Install MMseqs2 Version: f05f8c51d6e9c7c0b15fbd533e4b678303f50b3e
WORKDIR /tmp
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz &&\
    tar xvfz mmseqs-linux-sse41.tar.gz &&\
    rm mmseqs-linux-sse41.tar.gz &&\
    mv /tmp/mmseqs/bin/mmseqs /usr/local/bin &&\
    # remove useless files
    rm -r /tmp/mmseqs


# Install mafft 7.313
RUN rm /usr/bin/mafft  # remove mafft installed with bioperl
WORKDIR /tmp
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.313-with-extensions-src.tgz &&\
    tar xf mafft-7.313-with-extensions-src.tgz &&\
    rm mafft-7.313-with-extensions-src.tgz
WORKDIR /tmp/mafft-7.313-with-extensions/core
RUN make clean &&\
    make &&\
    make install


# Install FastTree version 2.1.11 Double precision (No SSE3)
WORKDIR /tmp
RUN wget http://www.microbesonline.org/fasttree/FastTree.c &&\
    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP FastTree.c -lm &&\
    ln -s /tmp/FastTreeMP /usr/local/bin


## Install FastME FastME 2.1.6.1
WORKDIR /tmp
RUN apt-get install -y automake  &&\
    git clone https://gite.lirmm.fr/atgc/FastME.git
WORKDIR /tmp/FastME
RUN ./configure &&\
    make &&\
    make install


# Install quicktree
WORKDIR /tmp
RUN git clone https://github.com/tseemann/quicktree
WORKDIR quicktree
RUN make &&\
    ln -s /tmp/quicktree/quicktree /usr/local/bin


# Install iqtree
WORKDIR /tmp
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz
RUN tar -xzf iqtree-1.6.12-Linux.tar.gz &&\
    rm iqtree-1.6.12-Linux.tar.gz &&\
    ln -s /tmp/iqtree-1.6.12-Linux/bin/iqtree /usr/local/bin


# Install iqtree2
WORKDIR /tmp
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.6/iqtree-2.0.6-Linux.tar.gz
RUN tar -xzf iqtree-2.0.6-Linux.tar.gz &&\
    rm iqtree-2.0.6-Linux.tar.gz &&\
    ln -s /tmp/iqtree-2.0.6-Linux/bin/iqtree2 /usr/local/bin


# Install PanACoTA
WORKDIR /tmp
RUN git clone https://gitlab.pasteur.fr/aperrin/pipeline_annotation.git
WORKDIR /tmp/pipeline_annotation
RUN git checkout dev &&\
    ./make


ENTRYPOINT ["/usr/local/bin/PanACoTA"]
CMD ['-h']
