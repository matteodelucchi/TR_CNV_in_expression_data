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
        -i1 | --input1 )          shift
                                infile1=$1
                                ;;
        -i2 | --input2 )          shift
                                infile2=$1
                                ;;
        -i3 | --input3 )          shift
                                infile3=$1
                                ;;
        -o | --output )         shift
                                outfile=$1
                                ;; 
    esac
    shift
done

echo "Processing:  ${infile1}, ${infile2} and ${infile3}"

# ----------------------------
# Short Tandem Repeat Calling
# ----------------------------
# https://hipstr-tool.github.io/HipSTR-tutorial/
# Run HipSTR
echo "Running HipSTR..."
/dataT/dlc/programs/HipSTR/HipSTR --bams       ${infile1},${infile2},${infile3} \
                                   --fasta     /dataT/dlc/data/hg38/all_chroms.fa \
                                   --regions   /dataT/dlc/data/hg38.hipstr_reference.bed \
                                   --str-vcf   ${outfile}.vcf.gz \
                                   --log       ${outfile}.log \
                                   --viz-out   ${outfile}.viz.gz \
                                   --min-reads 25 --def-stutter-model \

# ----------------------------
# Generate Consensus Sequence
# ----------------------------
# https://samtools.github.io/bcftools/howtos/consensus-sequence.html
echo "Generating consensus sequence..."
bcftools index ${outfile}.vcf.gz
cat /dataT/dlc/data/hg38/all_chroms.fa | bcftools consensus ${outfile}.vcf.gz > ${outfile}.consensus.fa

# ----------------------------
# Extract Coding Sequence
# ----------------------------
echo "WARNING: Extracting Coding Sequence is not yet implemented."
# Cufflinks/EMBOSS transeq
# extractfeat -type CDS \
#             -join \
#             -sequence test2.consensus.fa \
#             -outseq test2.features.fasta \
#             /dataT/dlc/data/hg38/all_chroms.fa
# extractfeat test2.consensus.fa test2.features.fasta
# NOTE: This doesn't work (both of them)!!!!! Keep on with all sequences, no only coding sequence
# Extract features from sequence(s)
# Warning: No sequences written to output file 'test2.features.fasta'

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
echo "Finished. Results in ${outfile}.blast.tab"

# # ----------------------------
# # Tandem Repeat Annotation (TRAL)
# # ----------------------------
# # TODO

