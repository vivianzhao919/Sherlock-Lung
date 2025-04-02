# RNA-seq data preprocess and quantification pipeline

This repository provides scripts for alignment and quantification of RNA-seq data for the Sherlock-Lung study (Genome Build: GRCh38).

## Requirements  
- Sentieon software (license required)
- Samtools, HTSeq
- Human reference genome: hg38 (GRCh38)
- Transcript annotation: GENCODE v35

## Workflow Overview
- Step 1. Alignment and quality assessment
  `sentieon_rna_star_list.sh` takes a list of input files and generate bash commands for alignment and quality assessment steps.
  ` Usage:
  ```sh ~/scripts/sentieon_rna_star_list.sh \
        -i <List_of_input_files> \
        -d <Output_directory> \
        -o <Output_bash_file \
        -n <nthreads>
  
- [Usage](#usage)  
- [Features](#features)  
- [Configuration](#configuration)  
- [Contributing](#contributing)  
- [License](#license)  
- [Acknowledgments](#acknowledgments)  

## Installation  
Provide step-by-step instructions on how to install and set up your project.  

```bash
# Example installation command
git clone https://github.com/username/repository.git
cd repository
pip install -r requirements.txt
