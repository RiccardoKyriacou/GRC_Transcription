# "Evidence for Transcription and Horizontal Gene Transfer in Dipteran Germline-Restricted Chromosomes" (2025)

### Overview

This directory contains all code for the paper 
"Evidence for Transcription and Horizontal Gene Transfer in Dipteran Germline-Restricted Chromosomes" (in press)
---

### Directory Summary

| Directory | Description |
|-----------|-------------|
| `01_RNAseq_Mapping` | Scripts to align RNA-seq reads to the reference _B. coprophila_ genome using STAR and calculate TPM values with StringTie. |
| `02_Intergenic_mapping` | Scripts to generate intergenic GFF files, map reads to intergenic regions, and calculate intergenic TPMs. |
| `03_Expressed_GRC_genes` | Scripts to identify expressed GRC-linked genes, BLAST GRC genes against core genome, and perform homology analysis against NCBI NR proteins and RepBase. |
| `04_HGT` | Scripts to investigate putative horizontal gene transfer (HGT) regions in GRC2. |
| `05_Urban_et_al_re-analysis` | Scripts to re-analyse pooled embryo RNA-seq data from Urban et al. (2021). |
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
