# 2. Mapping reads to intergenic regions 

First, we can use the following python script as so to generate intergenic GFF/GTF files, providing the minimum distance from genes and the minimum/maximum length for the intergenic regions 
```
python3 01_get_intergenic_GFF3.py -g ../../Annotations/bcop_core_GRC.gff3 --min_distance 500 --min_length 1000 --max_length 20000
```
Then we can run 02_intergenic_TPM.sh, running StringTie on the uniquley mapped BAM files from 01, this time using the intergenic gtf to calculate TPMs
```
for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename "$file" "_uniquely_mapped.bam")
    output_gtf="${base}_intergenic.gtf"
    echo "Calculating TPM for $file"
    stringtie "$file" -p 16 -o "$output_gtf" -G bcop_core_GRC.intergenic.gtf
done
```
Then we can run the next python script in the same directory where we ran StringTie to generate a combined TPM file 
```
pyhton 03_get_intergenic_TPM.py -t . 
```
Finally, the R script 04_intergenic_TPM_deconvolution.R visualises and performs statistics to generate a TPM cutoff value for expression, based on this intergenic mapping rate
