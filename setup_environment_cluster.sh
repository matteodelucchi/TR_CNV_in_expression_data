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
mkdir /dataT/dlc/data/hg19
cd /dataT/dlc/data/hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xzvf chromFa.tar.gz
cat chr*.fa > all_chroms.fa
samtools faidx all_chroms.fa

mkdir /dataT/dlc/data/hg38
cd /dataT/dlc/data/hg38
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz' -O chromFa.tar.gz
tar -xzvf chromFa.tar.gz
cat ./chroms/chr*.fa > ./all_chroms.fa
samtools faidx all_chroms.fa

# download regions .bed file with information about each STR.
# Here we use the prepared one from HipSTR.
# But this is a tutorial to create a new one: https://github.com/HipSTR-Tool/HipSTR-references/blob/master/human/human_reference.md
cd /dataT/dlc/data/
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh38.hipstr_reference.bed.gz
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz
gunzip *.hipstr_reference.bed.gz

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

# TODO: Add workaround for HMMER4-devel

# ----------------------------
# Install Eigen
# ----------------------------
cd /dataT/dlc/programs
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
tar xjf 3.3.7.tar.bz2
rm 3.3.7.tar.bz2 
cd eigen-eigen-323c052e1731/w
mkdir build_dir
cd build_dir
cmake ../

# ----------------------------
# Install Cufflinks
# ----------------------------
module load Boost/1.67.0-ibmatcuda-10.10.1
git clone https://github.com/cole-trapnell-lab/cufflinks.git
cd cufflinks
autoreconf --install
./configure --prefix=/dataT/dlc/programs/cufflinks \
            --with-boost=$MODULEPATH \
            --with-eigen=/dataT/dlc/programs/eigen-eigen-323c052e1731/build_dir \
            --with-bam=/dataT/dlc/samtools-1.9
# NOTE: FAILED INSTALLATION
# checking for bamlib... 
# configure: error: We could not detect the bam libraries (version
# or higher). If you have a staged bam library (still not installed) 
# please specify $BAM_ROOT in your environment and do not give a PATH 
# to --with-bam option.  If you are sure you have bam installed, 
# then check your version number looking in <bam/version.hpp>. 
# See http://randspringer.de/bam for more documentation.

# make
# make install
# # test with:
# wget http://cufflinks.cbcb.umd.edu/downloads/test_data.sam
# cufflinks ./test_data.sam

# wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
# tar zxvpf cufflinks-2.2.1.Linux_x86_64.tar.gz
# rm cufflinks-2.2.1.Linux_x86_64.tar.gz
# cd cufflinks-2.2.1.Linux_x86_64/
# echo -e "\nexport PATH=$PATH:/dataT/dlc/programs/cufflinks-2.2.1.Linux_x86_64
# " >> ~/.bashrc
# wget http://cufflinks.cbcb.umd.edu/downloads/test_data.sam
# cufflinks ./test_data.sam

# ----------------------------
# Install TransDecoder
# ----------------------------
cd /dataT/dlc/programs
wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
tar zxvpf TransDecoder-v5.5.0.tar.gz
rm TransDecoder-v5.5.0.tar.gz
cd TransDecoder-TransDecoder-v5.5.0/
autoconf
# NOTE: FAILED INSTALLATION requires Cufflinks (see above)
# ./configure
# # --build=ppc64le
# make clean
# make
# make install


# # other option with Conda
# module load Anaconda3/2018.12-rhel-7
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda install transdecoder
# conda update transdecoder
# NOTE: is not in the conda package library


# ----------------------------
# Install Picard
# ----------------------------
cd /dataT/dlc/programs
# git clone https://github.com/broadinstitute/picard.git
# cd picard/
# ./gradlew shadowJar
# java -jar build/libs/picard.jar
# ./gradlew clean
# ./gradlew test
# NOTE: didn't pass the tests

# mkdir picard
# cd ./picard/
# wget https://github.com/broadinstitute/picard/releases/download/2.20.2/picard.jar
# echo -e "\nexport PICARD=/dataT/dlc/programs/picard/picard.jar" >> ~/.bashrc
# NOTE: Doesn't work due to processor architecture incompatibilities...
# INFO	2019-06-27 10:59:29	SamToFastq	

# ********** NOTE: Picard's command line syntax is changing.
# **********
# ********** For more information, please see:
# ********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
# **********
# ********** The command line looks like this in the new syntax:
# **********
# **********    SamToFastq -I /dataT/dlc/tcga/colorectal/bam_files/primary_tumor/TCGA-A6-2685-01A-01D-1408-10_Illumina_gdc_realn.bam -FASTQ test_de_novo.fastq
# **********


# 10:59:30.037 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/ibm/gpfs-dataT/dlc/programs/picard/picard.jar!/com/intel/gkl/native/libgkl_compression.so
# 10:59:30.065 WARN  NativeLibraryLoader - Unable to load libgkl_compression.so from native/libgkl_compression.so (/tmp/dlc/libgkl_compression1486417447529522737.so: /tmp/dlc/libgkl_compression1486417447529522737.so: cannot open shared object file: No such file or directory (Possible cause: can't load AMD 64-bit .so on a Power PC 64 LE-bit platform))
# 10:59:30.066 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/ibm/gpfs-dataT/dlc/programs/picard/picard.jar!/com/intel/gkl/native/libgkl_compression.so
# 10:59:30.071 WARN  NativeLibraryLoader - Unable to load libgkl_compression.so from native/libgkl_compression.so (/tmp/dlc/libgkl_compression7237514213767117405.so: /tmp/dlc/libgkl_compression7237514213767117405.so: cannot open shared object file: No such file or directory (Possible cause: can't load AMD 64-bit .so on a Power PC 64 LE-bit platform))
# [Thu Jun 27 10:59:30 CEST 2019] SamToFastq INPUT=/dataT/dlc/tcga/colorectal/bam_files/primary_tumor/TCGA-A6-2685-01A-01D-1408-10_Illumina_gdc_realn.bam FASTQ=test_de_novo.fastq    OUTPUT_PER_RG=false COMPRESS_OUTPUTS_PER_RG=false RG_TAG=PU RE_REVERSE=true INTERLEAVE=false INCLUDE_NON_PF_READS=false CLIPPING_MIN_LENGTH=0 READ1_TRIM=0 READ2_TRIM=0 INCLUDE_NON_PRIMARY_ALIGNMENTS=false VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
# [Thu Jun 27 10:59:30 CEST 2019] Executing as dlc@zhcc068 on Linux 4.14.0-115.6.1.el7a.ppc64le ppc64le; OpenJDK 64-Bit Server VM 1.8.0_201-b09; Deflater: Jdk; Inflater: Jdk; Provider GCS is not available; Picard version: 2.20.2-SNAPSHOT
# 10:59:30.211 WARN  IntelDeflaterFactory - IntelInflater is not supported, using Java.util.zip.Inflater
# [Thu Jun 27 10:59:30 CEST 2019] picard.sam.SamToFastq done. Elapsed time: 0.01 minutes.
# Runtime.totalMemory()=81264640
# To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp
# Exception in thread "main" picard.PicardException: Input contains paired reads but no SECOND_END_FASTQ specified.
# 	at picard.sam.SamToFastq.handleRecord(SamToFastq.java:303)
# 	at picard.sam.SamToFastq.doWork(SamToFastq.java:192)
# 	at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:295)
# 	at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:103)
# 	at picard.cmdline.PicardCommandLine.main(PicardCommandLine.java:113)



# ----------------------------
# Install Bedtools
# ----------------------------
# cd /dataT/dlc/programs
# wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
# tar -zxvf bedtools-2.28.0.tar.gz
# cd bedtools2
# make
# echo -e "\nexport PATH=$PATH:/dataT/dlc/programs/bedtools2/bin" >> ~/.bashrc
# NOTE: Gives same error on multiple tested files:
# [dlc@zhcc054 .../dlc/de-novo_assembly]$ bamToFastq -i /dataT/dlc/tcga/colorectal/bam_files/primary_tumor/TCGA-A6-2671-01A-01D-1408-10_hg19_Illumina_gdc_realn.bam -fq test_de-novo.fq
# [E::bgzf_read] Read block operation failed with error 2 after 282 of 411 bytes

# ----------------------------
# Install Hydra
# ----------------------------
# cd /dataT/dlc/programs
# wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/hydra-sv/Hydra.v0.5.3.tar.gz
# tar -zxvf Hydra.v0.5.3.tar.gz 
# cd Hydra-Version-0.5.3/
# make clean
# make all
# echo -e "\nexport hydra=/dataT/dlc/programs/Hydra-Version-0.5.3/bin" >> ~/.bashrc
# NOTE: installation successfull but doesn't work though...