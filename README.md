# Tandem Repeat Copy Number Variation in \n Expression Data
This repository contains pipelines to preprocess RNAseq `.bam` files for the analysis with [TRAL](https://github.com/acg-team/tral).

## Setup and Installation of Dependencies
Download this repository
```bash
git clone https://github.com/matteodelucchi/TR_CNV_in_expression_data.git
```

Set up the directory struture for the files and install all necessary programs by sourcing the respective file for installation in a local environment or a high-performance cluster environment
```bash
chmod 779 setup_environment_cluster.sh
./setup_environment_cluster.sh
```

## Run Pipelines
The two proposed pipelines are provided in 
* `reference-based_mapping_pipeline.sh`
* `de-novo_assembly.sh`

They are on early reasearch stage - don't run them as executables.
Go through each line very carefully step by step.

