# RNA-seq data preprocess and quantification pipeline

This repository provides scripts for alignment and quantification of RNA-seq data for the Sherlock-Lung study (Genome Build: GRCh38).

## Requirements  
- Sentieon software (license required)
- Samtools, HTSeq
- Human reference genome: hg38 (GRCh38)
- Transcript annotation: GENCODE v35

## Workflow Overview
### Step 1. Alignment and quality assessment
`STAR_alignment.sh` aligns FASTQ files from RNA-seq data to genome files using STAR and performs quality assessment using Picard Tools.
### Usage:
  ```
  sh STAR_alignment.sh -s <Sample_name> -t <raw/trim/plain> -i <FASTQ_folder> -o <Output_folder> -n <nthreads>
  ```
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
  ```

### Step 2. Transcript quantification
`htseq_quantification.sh` takes BAM files aligned to the genome and generate
### Usage:
  ```
  sh htseq_quantification.sh sh -s <Sample_name> -i <BAM_folder> -o <Output_folder> -m <memory>
  ```
<Sample_name>: Input sample file name. The name of FASTQ file should contain Sample_name.  
<raw/trim/plain>: Type of FASTQ files. The file names of input FASTQ should match these patterns:
  - A "raw" file ends with "1/2.fastq.gz".
  - A "trim" file ends with "1/2*trimmed.fastq.gz".
  - A "plain" file ends with "1/2.fastq".
<BAM_folder>: Directory of the input BAM files.
<Output_folder>: Directory of the output quantification data.
<nthreads>: Memery required.

### Example:
  ```
  sh STAR_alignment.sh -s CSP107355 -t raw -i /data/Sherlock_Lung/BAM -o/data/Sherlock_Lung/htseq -m 16g
