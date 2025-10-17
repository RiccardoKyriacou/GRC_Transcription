# 3. Generating expressed GRC data 

For a GRC-linked gene to be considered expressed in this study, it must satisfy two conditions
- 1) Have an TPM > cutoff generated in /02_intergenic_TPM/ in at least 2 out of the three replicate libraries
  2) Not see any such expression in somatic libraries
  3) Not have a highly similar core-chromosome paralogue
 
Hence these scripts mainly sort the data to identify GRC-linked genes in this manner

First we run the script 01_get_expressed_genes.py. 
This script 
- 1) Generates a table of expressed genes only where expressed genes = genes with a TPM > intergenic_trheshold present in 2 or more samples
2) Parses gene IDs for expressed GRC genes (these are the genes we want to investigate)
3) Using gene ID's parses out the gene sequences from B_cor GRC genome
4) Outputs three FASTA files for GRC genes expressed in soma, germ, and both
5) Parses out all core genes from annotated B_cop genome
6) Produces combined fasta file contining all core chromosome genes
```
python 01_get_expressed_genes_for_BLAST.py -g ../Annotations/idBraCopr2.1.primary.masked -a ../Annotations/bcop_core_GRC.gff3 -o ../01_STAR_TPM/ -t 0.22 
```
