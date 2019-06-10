#!/bin/bash

# Script to setup directory structure and installation of all necessary programs
# necessary to run reference-based_mapping_pipeline.sh and de-novo_assembly_pipeline.sh

# ----------------------------
# Build Directory Structure
# ----------------------------
mkdir ./data_TCGA-COAD_refseq_bam # dbGAP Data from IBM
mkdir ./prog_HipSTR # install HipSTR with all it's dependencies here


# ----------------------------
# Install samtools
# ----------------------------
# https://samtools.github.io/bcftools/howtos/install.html

# ----------------------------
# Install HipSTR
# ----------------------------
# https://github.com/tfwillems/HipSTR#requirements
apt install make g++ zlib1g-dev libhts-dev libbz2-dev liblzma-dev

cd ./prog_HipSTR

git clone https://github.com/HipSTR-Tool/HipSTR
cd HipSTR
make

# ----------------------------
# Install EMBOSS
# ----------------------------
# use new, parallelizable version:
# https://github.com/feliixx/gotranseq

# ----------------------------
# Install BLAST+
# ----------------------------
# https://www.ncbi.nlm.nih.gov/books/NBK52640/
# or here specifically for redhat:
# https://www.ncbi.nlm.nih.gov/books/NBK279671/