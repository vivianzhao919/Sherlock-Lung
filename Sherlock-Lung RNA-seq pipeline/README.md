# RNA-seq data preprocess and quantification pipeline

This repository provides scripts for alignment and quantification of RNA-seq data for the Sherlock-Lung study (Genome Build: GRCh38).

## Requirements  
- Sentieon software (license required)
- Samtools, HTSeq
- Human reference genome: hg38 (GRCh38)
- Transcript annotation: GENCODE v35

## Workflow Overview
### Step 1. Alignment and quality assessment
`STAR_alignment.sh` takes a list of input files and generate bash commands for alignment and quality assessment steps.
### Usage:
  ```
  sh STAR_alignment.sh -s <Sample_name> -t <raw/trim/plain> -i <FASTQ_folder> -o <Output_folder> -n <nthreads>

<Sample_name>: Input sample file name. The name of FASTQ file should contain Sample_name.  
<raw/trim/plain>: Type of FASTQ files. The file names of input FASTQ should match these patterns:
  - A "raw" file ends with "1/2.fastq.gz".
  - A "trim" file ends with "1/2*trimmed.fastq.gz".
  - A "plain" file ends with "1/2.fastq".
<FASTQ_folder>: Directory of the input FASTQ files.
<Output_folder>: Directory of the output BAM and QC files.
<nthreads>: Number of CPU threads.

### Example:
  ```
  sh STAR_alignment.sh -s CSP107355 -t raw -i /data/Sherlock/FASTQ -o /data/Sherlock_Lung/BAM -n 4
  

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
