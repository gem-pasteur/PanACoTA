from ubuntu:22.04

# setting these ARGs as global variables, so they persist across all build layers (but not in final production image)
ARG PROKKA_VER="1.14.6"
ARG BEDTOOLS_VER="2.29.0"
ARG BARRNAP_VER="0.9"
ARG BLAST_VER="2.10.1"
ARG MASH_VER="v2.2"
ARG PRODIGAL_VER="2.6.3"
ARG MAFFT_VER="7.313"
ARG FASTME_VER="2.1.6.1"

# re-instantiating for the app build layer
ARG PROKKA_VER
ARG BEDTOOLS_VER
ARG BARRNAP_VER
ARG BLAST_VER
ARG MASH_VER
ARG PRODIGAL_VER
ARG MAFFT_VER
ARG FASTME_VER

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
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz &&\
    tar zxvpf ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz &&\
    cp /install_dir/ncbi-blast-${BLAST_VER}+/bin/makeblastdb /usr/local/bin/ &&\
    cp /install_dir/ncbi-blast-${BLAST_VER}+/bin/blastp /usr/local/bin/ &&\
    rm ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz

# Install mash
WORKDIR /install_dir
RUN wget https://github.com/marbl/Mash/releases/download/${MASH_VER}/mash-Linux64-${MASH_VER}.tar &&\
    tar -xf mash-Linux64-${MASH_VER}.tar &&\
    rm mash-Linux64-${MASH_VER}.tar &&\
    mv /install_dir/mash-Linux64-${MASH_VER}/mash /usr/local/bin &&\
    rm -r mash-Linux64-${MASH_VER}



# Install dependencies for prokka:
# install tbl2asn manually since the one included with prokka is expired. 
# Probably will have to do again in Dec 2021 unless Torsten removes it from prokka
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz && \
    gunzip linux64.tbl2asn.gz && \
    mv linux64.tbl2asn /usr/bin/tbl2asn && \
    chmod 755 /usr/bin/tbl2asn
# Install prokka
# RUN git clone https://github.com/tseemann/prokka.git
# RUN /install_dir/prokka/bin/prokka --setupdb &&\
#     ln -s /install_dir/prokka/bin/prokka /usr/local/bin
# # install tbl2asn (used by prokka)
# RUN wget -O tbl2asn.gz https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz &&\
#     gunzip tbl2asn.gz &&\
#     chmod +x tbl2asn &&\
#     ln -s /install_dir/tbl2asn /usr/local/bin
RUN apt-get install -y prokka

WORKDIR /install_dir
# RUN DEBIAN_FRONTEND="noninteractive" apt install -y\
#         libdatetime-perl \
#         libxml-simple-perl \
#         libdigest-md5-perl \
#         hmmer \
#         default-jre \
#         bioperl
# # Install hmmer
# RUN echo yes | cpan Bio::SearchIO::hmmer
# # Install bioperl
# RUN echo yes | cpan Bio::Perl
# # Install prokka
RUN wget https://github.com/tseemann/prokka/archive/refs/tags/v${PROKKA_VER}.tar.gz
RUN tar -xzf v${PROKKA_VER}.tar.gz &&\
    rm -rf v${PROKKA_VER}.tar.gz 
RUN /install_dir/prokka-${PROKKA_VER}/bin/prokka --setupdb &&\
    ln -s /install_dir/prokka-${PROKKA_VER}/bin/prokka /usr/local/bin

# # install tbl2asn (used by prokka)
# RUN wget -O tbl2asn.gz https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz &&\
#     gunzip tbl2asn.gz &&\
#     chmod +x tbl2asn &&\
#     ln -s /install_dir/tbl2asn /usr/local/bin
# more Examples/commands/2-Annotate.sh
# PanACoTA annotate -d Examples/genomes_init -r Examples/2-res-prokka -l Examples/input_files/list_genomes.lst -n GENO --l90 3 --nbcont 10
 

# Install prodigal
RUN apt-get install prodigal
# WORKDIR /install_dir
# RUN wget https://github.com/hyattpd/Prodigal/archive/v${PRODIGAL_VER}.tar.gz &&\
#     tar -xzf v${PRODIGAL_VER}.tar.gz &&\
#     rm v${PRODIGAL_VER}.tar.gz
# WORKDIR /install_dir/Prodigal-${PRODIGAL_VER}
# RUN make 
# RUN make install

# Install MMseqs2 Version: f05f8c51d6e9c7c0b15fbd533e4b678303f50b3e
WORKDIR /install_dir
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz &&\
    tar xvfz mmseqs-linux-sse41.tar.gz &&\
    rm mmseqs-linux-sse41.tar.gz &&\
    mv /install_dir/mmseqs/bin/mmseqs /usr/local/bin &&\
    # remove useless files
    rm -r /install_dir/mmseqs

# # Install mafft 7.313
# RUN rm /usr/bin/mafft  # remove mafft installed with bioperl
# WORKDIR /install_dir
# RUN wget https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VER}-with-extensions-src.tgz &&\
#     tar xf mafft-7.313-with-extensions-src.tgz &&\
#     rm mafft-7.313-with-extensions-src.tgz
# WORKDIR /install_dir/mafft-7.313-with-extensions/core
# RUN make clean &&\
#     make &&\
#     make install

# # Install FastTree version 2.1.11 Double precision (No SSE3)
WORKDIR /install_dir
RUN wget http://www.microbesonline.org/fasttree/FastTree.c &&\
    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP FastTree.c -lm &&\
    ln -s /install_dir/FastTreeMP /usr/local/bin

# Install FastME FastME 2.1.6.1
WORKDIR /install_dir
RUN apt-get install -y automake  &&\
    git clone https://gite.lirmm.fr/atgc/FastME.git 
WORKDIR /install_dir/FastME/tarball
RUN tar xzf fastme-2.1.6.4.tar.gz &&\
    rm fastme-2.1.6.4.tar.gz &&\
<<<<<<< HEAD
    ln -s /install_dir/FastME/tarball/fastme-2.1.6.2/binaries/fastme-2.1.6.4-linux64-omp /usr/local/bin/fastme
=======
    ln -s /install_dir/FastME/tarball/fastme-2.1.6.4/binaries/fastme-2.1.6.2-linux64-omp /usr/local/bin/fastme
>>>>>>> dev

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
RUN wget https://github.com/gem-pasteur/PanACoTA/archive/refs/tags/v1.4.0.tar.gz
RUN tar -xzf v1.4.0.tar.gz
WORKDIR /install-dir/PanACoTA-1.4.0
RUN ./make

ENTRYPOINT ["/usr/local/bin/PanACoTA"]
CMD ['-h']