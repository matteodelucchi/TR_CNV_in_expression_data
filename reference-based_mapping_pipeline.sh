#!/bin/bash

# processes RNAseq .bam files according the pipeline for TR detection
# with TRAL.
# 
# Dependencies:
#   - see setup_environment_local.sh or setup_environment_cluster.sh
# 
# Arguments:
#   --input RNAseq.bam file
#   --output filename for the outputfile
#
# Returns:
#   file specified with --output <output>.fasta

# ----------------------------
# Read Input Arguments
# ----------------------------
while [ "$1" != "" ]; do
    case $1 in
        -i | --input )          shift
                                infile=$1
                                ;;
        -o | --output )         shift
                                outfile=$1
                                ;; 
    esac
    shift
done

echo "Processing: [$infile]"

# ----------------------------
# Mapped Read Data Filtering
# ----------------------------
echo "Mapped Read Data Filtering..."
samtools view -q 10 -b infile > aligned_reads.q10.bam

# ----------------------------
# Short Tandem Repeat Calling
# ----------------------------
# https://hipstr-tool.github.io/HipSTR-tutorial/
# Prepare for HipSTR input
echo "STR Calling"

echo "Prepare for HipSTR input (sort and index)"
samtools sort -o my_sorted.bam aligned_reads.q10.bam
samtools index > outfile.bam

echo "Finished. Results in [$outfile]"

# uncomment the chunck below




# samtools index > my_sorted.bam


# # Download human hg19 reference files in fasta format
# mkdir fasta
# cd fasta
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# tar -xzvf chromFa.tar.gz
# cat chr*.fa > all_chroms.fa
# samtools faidx all_chroms.fa
# cd ../

# # Download human hg19 STR-reference file
# wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz
# gunzip hg19.hipstr_reference.bed.gz

# # Run HipSTR
# ./HipSTR/HipSTR --bams      bams/ERR194147.bam,bams/ERR194160.bam,bams/ERR194161.bam,bams/SRR826427.bam,bams/SRR826428.bam,bams/SRR826448.bam,bams/SRR826463.bam,bams/SRR826465.bam,bams/SRR826467.bam,bams/SRR826469.bam,bams/SRR826471.bam,bams/SRR826473.bam
#                 --fasta     fasta/all_chroms.fa
#                 --regions   regions.bed
#                 --str-vcf   trio.marshfield.vcf.gz
#                 --log       trio.marshfield.log
#                 --viz-out   trio.marshfield.viz.gz
#                 --min-reads 25 --def-stutter-model

# # ----------------------------
# # Generate Consensus Sequence
# # ----------------------------
# # https://samtools.github.io/bcftools/howtos/consensus-sequence.html
# cat reference.fa | bcftools consensus calls.vcf.gz > consensus.fa

# # ----------------------------
# # Extract Coding Sequence
# # ----------------------------
# # Cufflinks/EMBOSS transeq
# extractfeat -type CDS -join -seqid -coords -retainids -matchdescstart -seqfile
# $DESTINATION -force -o $OUT_FILE $REFERENCE

# # ----------------------------
# # Translation DNA -> AA
# # ----------------------------
# # Emboss transeq
# transeq -sequence <input file> -outseq <output file> -frame 6 -clean
# # or with Gotranseq
# cpus=$( ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w )
# ./gotranseq --sequence <input file> --outseq <output file> --frame 6 --numcpu $cpus --clean

# # ----------------------------
# # Transcript Annotation
# # ----------------------------
# BLASTX
# # ----------------------------
# # Tandem Repeat Annotation (TRAL)
# # ----------------------------