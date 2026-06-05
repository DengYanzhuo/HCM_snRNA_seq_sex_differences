# HCM_snRNA_seq_sex_differences
Single-nucleus RNA-seq analysis of sex-dimorphic cellular and molecular dysregulation in human hypertrophic cardiomyopathy (HCM).

This repository contains the complete analysis code for the manuscript:

**Single-nucleus profiling of human hypertrophic cardiomyopathy uncovers sex-dimorphic cellular and molecular dysregulation**  
Yanzhuo Deng, Zhuoran Liang, Yihao Liu, Baofa Sun  

## Abstract

Hypertrophic cardiomyopathy (HCM) is a prevalent cardiovascular disease characterized by ventricular hypertrophy that can progress to heart failure. Despite well-recognized sex disparities in HCM, the cellular and molecular mechanisms underlying these differences remain incompletely understood. Here we performed single-nucleus transcriptome analysis of left ventricular tissues from male and female patients with HCM and sex-matched non-failing controls. Female patients exhibited stronger upregulation of hypertrophic and fibrotic programs than male patients, indicative of more pronounced pathological remodeling. In contrast, male patients showed greater impairment of metabolic homeostasis and enhanced cellular stress responses, whereas females demonstrated more prominent immune activation. Subpopulation analyses of cardiomyocytes, fibroblasts and macrophages further identified disease-associated subpopulations and sex-specific key regulatory genes. Moreover, intercellular communication revealed enhanced remodeling-associated signaling networks in females. Collectively, our findings delineate the cellular basis of sex-specific differences in HCM and identify candidate molecular targets that may inform the development of sex-tailored therapeutic strategies for precision cardiovascular medicine.

## Data availability

The single-nucleus RNA sequencing (snRNA-seq) data used in this study are publicly available from the **Broad Institute Single Cell Portal** under project ID **SCP1303**:

- Direct link: [https://singlecell.broadinstitute.org/single_cell/study/SCP1303/](https://singlecell.broadinstitute.org/single_cell/study/SCP1303/)

All data were downloaded and processed as described in the manuscript.

## System requirements and dependencies

All analyses were performed in **R (version 4.2.0 or higher)**. The following R packages are required:

| Package | Version |
|---------|---------|
| Seurat | 4.3.0 |
| clusterProfiler | 4.8.3 |
| GSVA | 2.0.7 |
| GSEA (clusterProfiler) | 4.8.3 |
| scMetabolism | 0.2.1 |
| AUCell | 1.22.0 |
| WGCNA | 1.72-1 |
| Monocle | 2.28.0 |
| CellChat | 1.6.1 |
| ggplot2 | 3.4.4 |
| ComplexHeatmap | 2.16.0 |
| UpSetR | 1.4.0 |
| limma | 3.62.2 |
