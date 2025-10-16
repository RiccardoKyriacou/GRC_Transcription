# GitHub Repository for "Evidence for Transcription and Horizontal Gene Transfer in Dipteran Germline-Restricted Chromosomes" (2025)

### Overview

#### This directory contains all code for the paper 
#### "Evidence for Transcription and Horizontal Gene Transfer in Dipteran Germline-Restricted Chromosomes" (in press)
---

### Directory Summary

1. **01_RNAseq_Mapping** Scripts to align RNA-seq reads to the reference _B. coprophila_ genome using the STAR genome aligner and calculate Transcript-per-million (TPM) values using StringTie.  
2. **02_Intergenic_mapping** Scripts to generate intergenic GFF, map reads to integrenic regions and caculate intergenic TPMs
3. **03_Expressed_GRC_genes** Scripts to get expressed GRC-linked genes, BLAST GRC genes agaisnt core genes/genome, and BLAST these genes against NCBI non-redundant proteins database and repbase to perform homology analysis.  
4. **04_HGT** Script used to investigate putatuve HGT region in GRC2
5. 4. **04_HGT** Script used to re-analyse pooled ambryo data from Urban et al (2021) 
---

### Directory Structure
```text
GRC_transcription/
├── 01_RNAseq_Mapping/
│   ├── 01_STAR_TPM_Nmax2
│   ├── 02_get_TPM_values.py
│   └── outputs/
│       ├── combined_TPM_only.tsv
│       ├── gene_overlap_0-4h_TPM.tsv
│       ├── gene_overlap_4-8h_TPM.tsv
│       ├── gene_overlap_adult_TPM.tsv
│       └── gene_overlap_late-larva-early-pupa_TPM.tsv
├── 02_Intergenic_mapping/
├── 03_Expressed_GRC_genes/
├── 04_GRC2_HGT/
├── 05_Urban_et_al_re-analysis/
├── figures/
└── README.md
```
