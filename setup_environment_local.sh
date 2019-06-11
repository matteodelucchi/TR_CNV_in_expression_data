#!/bin/bash

# Script to setup directory structure and installation of all necessary programs
# necessary to run reference-based_mapping_pipeline.sh and de-novo_assembly_pipeline.sh
# 
# Tested on:
# Red Hat Enterprise Linux Workstation release 7.6 (Maipo)
# Architecture: x86_64

# ----------------------------
# Build Directory Structure
# ----------------------------
mkdir ./data
mkdir ./data/TCGA-COAD_refseq_bam # dbGAP Data from IBM
mkdir ./programs
mkdir ./programs/HipSTR # install HipSTR with all it's dependencies here


# ----------------------------
# Install samtools
# ----------------------------
# http://www.htslib.org/download/
cd ~/TR_CNV_in_expression_data/programs/
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
sudo tar xjf samtools-1.9.tar.bz2 
rm samtools-1.9.tar.bz2 
cd samtools-1.9/
sudo ./configure
sudo make
sudo make install
samtools --version

cd ~/TR_CNV_in_expression_data/programs/
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
sudo tar xjf bcftools-1.9.tar.bz2 
rm bcftools-1.9.tar.bz2
cd bcftools-1.9/
sudo ./configure
sudo make
sudo make install
bcftools --version

cd ~/TR_CNV_in_expression_data/programs/
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
sudo tar xjf htslib-1.9.tar.bz2
rm htslib-1.9.tar.bz2
cd htslib-1.9/
sudo ./configure
sudo make
sudo make install

# ----------------------------
# Install HipSTR
# ----------------------------
# https://github.com/tfwillems/HipSTR#requirements
#apt install make g++ zlib1g-dev libhts-dev libbz2-dev liblzma-dev
sudo yum install zlib-devel libbz2-devel liblzma-devel
cd ~/TR_CNV_in_expression_data/programs/

git clone https://github.com/HipSTR-Tool/HipSTR
cd HipSTR
make
./HipSTR --help

# ----------------------------
# Install EMBOSS
# ----------------------------
# use new, parallelizable version:
# https://github.com/feliixx/gotranseq
# install go
cd ~/TR_CNV_in_expression_data/programs
wget https://dl.google.com/go/go1.12.5.linux-amd64.tar.gz
tar -xzf go1.12.5.linux-amd64.tar.gz 
sudo mv go /usr/local
sudo rm go1.12.5.linux-amd64.tar.gz 
export GOROOT=/usr/local/go
export GOPATH=$HOME/TR_CNV_in_expression_data
export PATH=$GOPATH/bin:$GOROOT/bin:$PATH
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
cd ~/TR_CNV_in_expression_data/programs
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar zxvpf ncbi-blast-2.9.0+-x64-linux.tar.gz
sudo rm -r ncbi-blast-2.9.0+-x64-linux.tar.gz
export PATH=$PATH:$HOME/ncbi-blast-2.9.0+/bin

# setup BLAST database
mkdir $HOME/blastdb
$export BLASTDB=$HOME/blastdb
cd $HOME/blastdb
ftp ftp.ncbi.nlm.nih.gov # Name: anonymous, Password: <mailadress>
cd blast/db
bin
get swissprot.tar.gz
bye
tar zxvpf swissprot.tar.gz
rm swissprot.tar.gz