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

echo "Processing:  ${infile}"

# ----------------------------
# Mapped Read Data Filtering
# ----------------------------
echo "Mapped Read Data Filtering..."
samtools view -q 10 -b ${infile} > ${outfile}.aligned_reads.q10.bam
samtools quickcheck ${outfile}.aligned_reads.q10.bam

# ----------------------------
# Short Tandem Repeat Calling
# ----------------------------
# https://hipstr-tool.github.io/HipSTR-tutorial/
# Prepare for HipSTR input
echo "STR Calling..."

echo "Prepare for HipSTR input (sort and index)..."
samtools sort -o ${outfile}.my_sorted.bam ${outfile}.aligned_reads.q10.bam
samtools quickcheck ${outfile}.my_sorted.bam

samtools index ${outfile}.my_sorted.bam
samtools quickcheck ${outfile}.my_sorted.bam

