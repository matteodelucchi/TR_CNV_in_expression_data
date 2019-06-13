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
samtools view -b -q 10 -b ${infile} > aligned_reads.q10.bam

# ----------------------------
# Short Tandem Repeat Calling
# ----------------------------
# https://hipstr-tool.github.io/HipSTR-tutorial/
# Prepare for HipSTR input
echo "STR Calling..."

echo "Prepare for HipSTR input (sort and index)..."
samtools sort -o my_sorted.bam aligned_reads.q10.bam
samtools index my_sorted.bam

# Run HipSTR
echo "Running HipSTR..."
/dataT/dlc/programs/HipSTR/HipSTR --bams      my_sorted.bam \
                                   --fasta     /dataT/dlc/data/all_chroms.fa \
                                   --regions   /dataT/dlc/data/hg19.hipstr_reference.bed \
                                   --str-vcf   ${outfile}.vcf.gz \
                                   --log       ${outfile}.log \
                                   --viz-out   ${outfile}.viz.gz \
                                   --min-reads 25 --def-stutter-model \

echo "Finished. Results in ${outfile}.vcf.gz"

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