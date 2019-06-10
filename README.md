# Tandem Repeat Copy Number Variation in Expression Data
This repository contains analysis pipelines to preprocess RNAseq `.bam` files for the analysis with [TRAL](https://github.com/acg-team/tral).

## Setup and Installation of Dependencies
Download this repository
```bash
mkdir TR_CNV_IN_EXPRESSION_DATA
cd TR_CNV_IN_EXPRESSION_DATA
git clone https://github.com/matteodelucchi/TR_CNV_IN_EXPRESSION_DATA
```

Set up the directory struture for the files and install all necessary programs by sourcing this file
```bash
chmod 779 setup_environment.sh
./setup_environment.sh
```

## Run Pipelines
The two proposed pipelines are provided in 
* `reference-base_mapping_pipeline.sh`
* `de-novo_assembly.sh`

They are on early reasearch stage - don't run them as executables.
Go through each line very carefully step by step.

