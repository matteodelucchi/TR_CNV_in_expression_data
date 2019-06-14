#!/bin/bash

# Script to setup directory structure and installation of all necessary programs
# necessary to run reference-based_mapping_pipeline.sh and de-novo_assembly_pipeline.sh
# 
# Tested on:
# Red Hat Enterprise Linux Server release 7.6 (Maipo)
# Architecture: ppc64le

# ----------------------------
# Build Directory Structure
# ----------------------------
mkdir /dataT/dlc/programs
mkdir /dataT/dlc/data/

# ----------------------------
# Install samtools
# ----------------------------
# http://www.htslib.org/download/
cd /dataT/dlc/programs/
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xjf samtools-1.9.tar.bz2 
rm samtools-1.9.tar.bz2 
cd samtools-1.9/

# vim ~/.bashrc # append: export PATH:/dataT/dlc/programs/bin"
# source ~/.bashrc 
echo -e "\nexport PATH:/dataT/dlc/programs/bin" >> ~/.bashrc
source ~/.bashrc

./configure --prefix=/dataT/dlc/programs
make
make install
samtools --version

cd /dataT/dlc/programs/
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar xjf bcftools-1.9.tar.bz2 
rm bcftools-1.9.tar.bz2
cd bcftools-1.9/
./configure --prefix=/dataT/dlc/programs
make
make install
bcftools --version

cd /dataT/dlc/programs/
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar xjf htslib-1.9.tar.bz2
rm htslib-1.9.tar.bz2
cd htslib-1.9/
./configure --prefix=/dataT/dlc/programs
make
make install

# ----------------------------
# Install HipSTR
# ----------------------------
# https://github.com/tfwillems/HipSTR#requirements
#apt install make g++ zlib1g-dev libhts-dev libbz2-dev liblzma-dev
#yum install zlib-devel libbz2-devel liblzma-devel

cd /dataT/dlc/programs/
git clone https://github.com/HipSTR-Tool/HipSTR
cd HipSTR
make
./HipSTR --help

# download reference proteome
cd /dataT/dlc/data/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xzvf chromFa.tar.gz
cat chr*.fa > all_chroms.fa
samtools faidx all_chroms.fa

# download regions .bed file with information about each STR.
# Here we use the prepared one from HipSTR.
# But this is a tutorial to create a new one: https://github.com/HipSTR-Tool/HipSTR-references/blob/master/human/human_reference.md
cd /dataT/dlc/data/
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz
gunzip -k hg19.hipstr_reference.bed.gz

echo -e "\nexport PATH=$PATH:/dataT/dlc/data/" >> ~/.bashrc
source ~/.bashrc

# ----------------------------
# Install EMBOSS
# ----------------------------
# EMBOSS for extractfeat
cd /dataT/dlc/programs
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
gunzip EMBOSS-6.6.0.tar.gz 
tar xvf EMBOSS-6.6.0.tar
rm -r EMBOSS-6.6.0.tar
cd EMBOSS-6.6.0/
./configure --build=ppc64le --prefix=/dataT/dlc/programs
make
make install

# use new, parallelizable version:
# https://github.com/feliixx/gotranseq
# install go
cd /dataT/dlc/programs
wget https://dl.google.com/go/go1.12.5.linux-amd64.tar.gz
tar -xzf go1.12.5.linux-amd64.tar.gz 
mv go /usr/local
rm -r go1.12.5.linux-amd64.tar.gz 
# export GOROOT=/usr/local/go
# export GOPATH=$HOME/TR_CNV_in_expression_data
# export PATH=$GOPATH/bin:$GOROOT/bin:$PATH
echo -e "\nexport GOROOT=/usr/local/go" >> ~/.bashrc
source ~/.bashrc
echo -e "\nexport GOPATH=$HOME/TR_CNV_in_expression_data" >> ~/.bashrc
source ~/.bashrc
echo -e "\nexport PATH=$GOPATH/bin:$GOROOT/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
go version
go env

# install gotranseq
go get -u "github.com/feliixx/gotranseq"
gotranseq --help

# ----------------------------
# Install BLAST+
# ----------------------------
# https://www.ncbi.nlm.nih.gov/books/NBK52640/
# or here specifically for redhat:
# https://www.ncbi.nlm.nih.gov/books/NBK279671/
# for cluster:
# https://github.com/ppc64le/build-scripts/blob/master/ncbi-blast/ncbi-blast_rhel_7.3.sh
cd /dataT/dlc/programs
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-src.tar.gz
tar zxvpf ncbi-blast-2.9.0+-src.tar.gz
rm ncbi-blast-2.9.0+-src.tar.gz
cd ncbi-blast-2.9.0+-src/c++
./configure --build=ppc64le #check the architecture with lscpu
cd ReleaseMT/build
make all_r

# export PATH=$PATH:/dataT/dlc/programs/ncbi-blast-2.9.0+/c++/ReleaseMT/bin
echo -e "\nexport PATH=$PATH:/dataT/dlc/programs/ncbi-blast-2.9.0+/c++/ReleaseMT/bin" >> ~/.bashrc
source ~/.bashrc

# setup BLAST database
mkdir /dataT/dlc/programs/blastdb
# export BLASTDB=/dataT/dlc/programs/blastdb
echo -e "\nexport BLASTDB=/dataT/dlc/programs/blastdb" >> ~/.bashrc
source ~/.bashrc

cd /dataT/dlc/programs/blastdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar zxvpf swissprot.tar.gz
rm swissprot.tar.gz

# ----------------------------
# Install TRAL
# ----------------------------
cd /dataT/dlc/programs
mkdir tral_ext_software
mkdir tral_repository
cd /dataT/dlc/programs/tral_repository
git clone https://github.com/acg-team/tral.git

module load python/3.5-ibmatcuda-10.10.1 
python3 -m venv env_tral
source env_tral/bin/activate
cd /dataT/dlc/programs/tral_repository/tral/easy_setup
./setupTRAL.sh setup # yes to p-value download
