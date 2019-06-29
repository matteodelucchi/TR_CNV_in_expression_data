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
# Convert BAM to FASTQ
# ----------------------------
echo "Converting BAM to FASTQ..."
# With Picard (see installation issues due to pcc64le architecture)
# java -jar $PICARD SamToFastq \
#      I=${infile} \
#      FASTQ=${outfile}.fastq
# NOTE: eventhough error messages due to pcc64le are raised, a fastq file seems to be created...

# With bamtools (see installation script: issues with reading file)
# bamToFastq -i ${infile} \
#            -fq ${outfile}.fastq
# NOTE: eventhough error message due to reading file, a .fq file is produced.
# Comparing the sizes of the one from picard and the one from bamtools, results in a filesize of 0.2MB for picard and 25GB for bedtools which we assume to be more reasonable
# Proceeding with bedtools

# With hydra
# /dataT/dlc/programs/Hydra-Version-0.5.3/bin/bamToFastq -i ${infile} \
#                                                        -fq1 ${outfile}.fastq
# requires two fastq files. Probably for both consensus strands? 
# Didn't work though...