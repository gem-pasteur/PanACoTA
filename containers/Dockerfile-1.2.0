from ubuntu:20.04

# Update apt-get packages
RUN apt-get update &&\
    apt-get -y upgrade


# Install packages needed and update pip
RUN apt-get install -y \
        wget \
        python3 \
        python3-pip \
        git
# Upgrade pip
RUN pip3 install --upgrade pip
RUN mkdir /install_dir


# Update makeblastdb and blastp for prokka
WORKDIR /install_dir
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz &&\
    tar zxvpf ncbi-blast-2.10.1+-x64-linux.tar.gz &&\
    cp /install_dir/ncbi-blast-2.10.1+/bin/makeblastdb /usr/local/bin/ &&\
    cp /install_dir/ncbi-blast-2.10.1+/bin/blastp /usr/local/bin/ &&\
    rm ncbi-blast-2.10.1+-x64-linux.tar.gz


# Install mash
WORKDIR /install_dir
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar &&\
    tar -xf mash-Linux64-v2.2.tar &&\
    rm mash-Linux64-v2.2.tar &&\
    mv /install_dir/mash-Linux64-v2.2/mash /usr/local/bin &&\
    rm -r mash-Linux64-v2.2


# Install barrnap
WORKDIR /install_dir
RUN wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz &&\
    tar -xf 0.8.tar.gz &&\
    rm 0.8.tar.gz &&\
    mv /install_dir/barrnap-0.8/bin/barrnap /usr/local/bin  &&\
    # Remove heavy useless files
    rm -r /install_dir/barrnap-0.8/examples /install_dir/barrnap-0.8/build/*.aln


# Install prodigal
WORKDIR /install_dir
RUN wget https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz &&\
    tar -xzf v2.6.3.tar.gz &&\
    rm v2.6.3.tar.gz
WORKDIR /install_dir/Prodigal-2.6.3
RUN make &&\
    make install


# Install dependencies for prokka:
WORKDIR /install_dir
RUN DEBIAN_FRONTEND="noninteractive" apt install -y\
        libdatetime-perl \
        libxml-simple-perl \
        libdigest-md5-perl \
        hmmer \
        default-jre \
        bioperl
# Install hmmer
RUN echo yes | cpan Bio::SearchIO::hmmer
# Install bioperl
RUN echo yes | cpan Bio::Perl

# Install prokka
RUN git clone https://github.com/tseemann/prokka.git
RUN /install_dir/prokka/bin/prokka --setupdb &&\
    ln -s /install_dir/prokka/bin/prokka /usr/local/bin
# install tbl2asn (used by prokka)
RUN wget -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz &&\
    gunzip tbl2asn.gz &&\
    chmod +x tbl2asn &&\
    ln -s /install_dir/tbl2asn /usr/local/bin


# Install MMseqs2 Version: f05f8c51d6e9c7c0b15fbd533e4b678303f50b3e
WORKDIR /install_dir
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz &&\
    tar xvfz mmseqs-linux-sse41.tar.gz &&\
    rm mmseqs-linux-sse41.tar.gz &&\
    mv /install_dir/mmseqs/bin/mmseqs /usr/local/bin &&\
    # remove useless files
    rm -r /install_dir/mmseqs


# Install mafft 7.313
RUN rm /usr/bin/mafft  # remove mafft installed with bioperl
WORKDIR /install_dir
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.313-with-extensions-src.tgz &&\
    tar xf mafft-7.313-with-extensions-src.tgz &&\
    rm mafft-7.313-with-extensions-src.tgz
WORKDIR /install_dir/mafft-7.313-with-extensions/core
RUN make clean &&\
    make &&\
    make install


# Install FastTree version 2.1.11 Double precision (No SSE3)
WORKDIR /install_dir
RUN wget http://www.microbesonline.org/fasttree/FastTree.c &&\
    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP FastTree.c -lm &&\
    ln -s /install_dir/FastTreeMP /usr/local/bin


# Install FastME FastME 2.1.6.1
WORKDIR /install_dir
RUN apt-get install -y automake  &&\
    git clone https://gite.lirmm.fr/atgc/FastME.git
WORKDIR /install_dir/FastME/tarball
RUN tar xzf fastme-2.1.6.2.tar.gz &&\
    rm fastme-2.1.6.2.tar.gz &&\
    ln -s /install_dir/FastME/tarball/fastme-2.1.6.2/binaries/fastme-2.1.6.2-linux64-omp /usr/local/bin/fastme

# Install quicktree
WORKDIR /install_dir
RUN git clone https://github.com/tseemann/quicktree
WORKDIR /install_dir/quicktree
RUN make &&\
    ln -s /install_dir/quicktree/quicktree /usr/local/bin


# Install iqtree
WORKDIR /install_dir
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz
RUN tar -xzf iqtree-1.6.12-Linux.tar.gz &&\
    rm iqtree-1.6.12-Linux.tar.gz &&\
    ln -s /install_dir/iqtree-1.6.12-Linux/bin/iqtree /usr/local/bin


# Install iqtree2
WORKDIR /install_dir
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.6/iqtree-2.0.6-Linux.tar.gz
RUN tar -xzf iqtree-2.0.6-Linux.tar.gz &&\
    rm iqtree-2.0.6-Linux.tar.gz &&\
    ln -s /install_dir/iqtree-2.0.6-Linux/bin/iqtree2 /usr/local/bin


# Install PanACoTA
WORKDIR /install-dir
RUN git clone https://gitlab.pasteur.fr/aperrin/pipeline_annotation.git
WORKDIR /install-dir/pipeline_annotation
RUN git checkout master
RUN ./make


ENTRYPOINT ["/usr/local/bin/PanACoTA"]
CMD ['-h']
