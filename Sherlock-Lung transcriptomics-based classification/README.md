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
### Sherlock_RNA_clinical_data.txt
Characteristics and clinical data of patients in the Sherlock cohort.
### GIS_clinical_data.txt
Characteristics and clinical data of patients in the GIS cohort.
### Sherlock_WGS_driver_mut_fusion_bin.txt
Driver mutations and fusions in the Sherlock cohort.
### GIS_driver_mut_fusion_bin.txt
Driver mutations and fusions in the GIS cohort.
### Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt
Predicted NS-LUAD expression subtypes and distance to the centroids of subtypes in the Sherlock cohort.
### GIS_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt
Predicted NS-LUAD expression subtypes and distance to the centroids of subtypes in the GIS cohort.
### Sherlock_RNA_TN_GIS_info.txt
Sample list of the Sherlock and GIS cohorts.

