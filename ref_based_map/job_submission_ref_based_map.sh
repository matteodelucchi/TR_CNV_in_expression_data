#!/bin/bash

# ----------------------------
# Get all specimen names
# ----------------------------
# specimens={find /dataT/dlc/tcga/colorectal/bam_files/blood_derived_normal/.  -printf "%f\n" | cut -c1-12dl}
# for x in /dataT/dlc/tcga/colorectal/bam_files/blood_derived_normal/*.bam; do
#   mkdir "${x%.*}" && mv "$x" "${x%.*}"
# done

# cd /dataT/dlc/tcga/colorectal/bam_files/blood_derived_normal
# for filename in ./*; do
#     echo "${filename}" | cut -c3-14;
# done   


cd /dataT/dlc/tcga/colorectal/bam_files/blood_derived_normal
for filename in ./*; do
    spec=$(echo "${filename##*/}" | cut -c1-12) 

    # ----------------------------
    # Create a directory for each specimen
    # ----------------------------
    echo '------------------------------------------------'
    echo 'create directory:' $spec
    if [ ! -d "/dataT/dlc/ref_based_map/$spec" ]; then
        mkdir "/dataT/dlc/ref_based_map/$spec"
    else
        echo "directoy already exist. Nothing to do."
    fi

    # ----------------------------
    # write .lsf file for each specimen
    # ----------------------------
    echo 'create ' $spec'.lsf'
    # Copy template file in specific directory
    cp /dataT/dlc/template.lsf /dataT/dlc/ref_based_map/$spec/$spec.lsf # TODO: not sure if Job array works well with this naming.... evtl. change to spec1, spec2,...

    # Insert BSUB commands in the header
    sed -e "9i#BSUB -e ""$spec".stderr.%J"" -i /dataT/dlc/ref_based_map/$spec/$spec.lsf
    sed -e "10i#BSUB -o ""$spec".stdout.%J"" -i /dataT/dlc/ref_based_map/$spec/$spec.lsf
    sed -e "11i#BSUB -J ""$spec""" -i /dataT/dlc/ref_based_map/$spec/$spec.lsf

    # Append lines specific to each specimen
    echo "./reference-based_map_sort_filter.sh --input /dataT/dlc/tcga/colorectal/bam_files/primary_tumor/$spec*.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                     --output $spec.1" >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "./reference-based_map_sort_filter.sh --input /dataT/dlc/tcga/colorectal/bam_files/blood_derived_normal/$spec*.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                     --output $spec.2" >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "./reference-based_map_sort_filter.sh --input /dataT/dlc/tcga/colorectal/bam_files/solid_tissue_normal/$spec*.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                     --output $spec.3" >> /dataT/dlc/ref_based_map/$spec/$spec.lsf

    echo "./reference-based_STR_calling_TR-annotation.sh \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                            --input1 $spec.1.my_sorted.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                            --input2 $spec.2.my_sorted.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                            --input3 $spec.3.my_sorted.bam \ " >> /dataT/dlc/ref_based_map/$spec/$spec.lsf
    echo "                                            --output $spec" >> /dataT/dlc/ref_based_map/$spec/$spec.lsf

    echo "Batch Job submission prepared."

    # ----------------------------
    # Append job to job submission file
    # ----------------------------
    echo 'add job to job submission file'
    queue="prod.long" # dev (2h), prod.short (1h), prod.med (12h) or prod.long (2d)
    echo "-q $queue /dataT/dlc/ref_based_map/$spec/$spec.lsf" >> /dataT/dlc/ref_based_map/job_submission_file.txt

done   

# ----------------------------
# Clean up bug (first filename which is processed is *)
# ----------------------------
rm -r /dataT/dlc/ref_based_map/'*'
sed '1d' /dataT/dlc/ref_based_map/job_submission_file.txt > tmpfile; mv tmpfile /dataT/dlc/ref_based_map/job_submission_file.txt # POSIX

# ----------------------------
# Submit Job Pack
# ----------------------------
echo '------------------------------------------------'
echo "Undo the mis-en-place with:     rm -r TCGA* job_submission_file.txt"
echo 'Want to submit the job pack now?'
echo 'You could do it later manually by:     bsub -pack /dataT/dlc/ref_based_map/job_submission_file.txt'
read -p "(y/n)" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    bsub -pack /dataT/dlc/ref_based_map/job_submission_file.txt
fi