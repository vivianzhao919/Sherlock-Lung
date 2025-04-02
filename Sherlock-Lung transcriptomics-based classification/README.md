# Sherlock-Lung transcriptomics-based classification

## Overview  
This repository provides scripts for the analysis of the RNA-seq data from the Sherlock lung study, including cell type deconvolution, transcriptomics-based classification using the 60-gene signature and prognostic analysis.

## Requirements  
- Sentieon software (license required)
- Samtools, HTSeq
- Human reference genome: hg38 (GRCh38)
- Transcript annotation: GENCODE v35
---
## Scripts
### Cell_type_deconvolution.R
Performs deconvolution of major cell types using CAMTHC and references from Maynard et. al. (PMID:32822576)  and deconvolution of stromal cells using Bisque and references from Lambrechts et. al. (PMID:29988129).
### Clanc_classification_CV.R: 
Estimates the minimum gene sets to recapitulate the NS-LUAD expression subtypes and the error rate of predicting NS-LUAD expression subtypes using the 60-gene signature in the Sherlock cohort.
### GIS_nomalization.R: 
Normalizes the RNA-seq data of the Sherlock (training) and GIS (testing) sets to perform classification.
### Clanc_classification.R: 
Performs classification of RNA-seq data using the 60-gene predictor in the Sherlock (training) and GIS (testing) sets.
### Multivariate_model_for_survival_prediction.R: 
Tests the prognostic prediction in the Sherlock and GIS data set using different sets of molecular and clinical features

---
## Data
