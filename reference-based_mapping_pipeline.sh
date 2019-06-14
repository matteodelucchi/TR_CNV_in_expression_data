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

# Run HipSTR
echo "Running HipSTR..."
/dataT/dlc/programs/HipSTR/HipSTR --bams      ${outfile}.my_sorted.bam \
                                   --fasta     /dataT/dlc/data/all_chroms.fa \
                                   --regions   /dataT/dlc/data/hg19.hipstr_reference.bed \
                                   --str-vcf   ${outfile}.vcf.gz \
                                   --log       ${outfile}.log \
                                   --viz-out   ${outfile}.viz.gz \
                                   --min-reads 25 --def-stutter-model \


# ----------------------------
# Generate Consensus Sequence
# ----------------------------
# https://samtools.github.io/bcftools/howtos/consensus-sequence.html
echo "Generating consensus sequence..."
cat /dataT/dlc/data/all_chroms.fa | bcftools consensus ${outfile}.vcf.gz > ${outfile}.consensus.fa

# ----------------------------
# Extract Coding Sequence
# ----------------------------
# Cufflinks/EMBOSS transeq
# extractfeat -type CDS \
#             -join \
#             test1.features.fasta \
#             /dataT/dlc/data/all_chroms.fa
# extractfeat test1.consensus.fa test1.features.fasta
# NOTE: This doesn't work!!!!! Keep on with all sequences, no only coding sequence

# # ----------------------------
# # Translation DNA -> AA
# # ----------------------------
# # Emboss transeq
echo "Translating DNA to AA..."
transeq -sequence ${outfile}.consensus.fa -outseq ${outfile}.aa.fasta -frame 6 -clean
# # or with Gotranseq
# cpus=$( ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w )
# gotranseq --sequence <input file> --outseq <output file> --frame 6 --numcpu $cpus --clean
# gotranseq --sequence test1.consensus.fa --outseq test1.aa.fasta --frame 6 --numcpu 4 --clean
# NOTE: gotranseq shows strange behaviour: prints the DNA seq on the screen... if put in a file ...>out.gotranseq it takes ages and if interupted its a binary file...


# # ----------------------------
# # Transcript Annotation
# # ----------------------------
echo "Transcript annotation with BLASTp..."
/dataT/dlc/programs/ncbi-blast-2.9.0+-src/c++/ReleaseMT/bin/blastp \
    -query ${outfile}.aa.fasta \
    -db swissprot.00 \
    -out ${outfile}.blast.tab \
    -outfmt 7

# ----------------------------
# Tandem Repeat Annotation (TRAL)
# ----------------------------
# TODO

echo "Finished. Results in ${outfile}.blast.tab"