#!/bin/bash -l
#SBATCH --job-name=download
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --export=ALL
#SBATCH --time=00-120:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=16gb
#SBATCH --output=download.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

###
# NOTE: All RNA-seq libraries from:
# Urban, J.M., Foulk, M.S., Bliss, J.E. et al. 
# High contiguity de novo genome assembly and DNA modification analyses for the fungus fly, Sciara coprophila, using single-molecule sequencing. 
# BMC Genomics 22, 643 (2021). https://doi.org/10.1186/s12864-021-07926-2
###

set -euo pipefail

srx_list=(SRX6716708 SRX6716709 SRX6716710 SRX6716711 SRX6716704 SRX6716705)

outdir="embryo_fastqs"
mkdir -p "$outdir"

for srx in "${srx_list[@]}"; do
  echo "Fetching runs for $srx ..."
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srx}&result=read_run&fields=run_accession,fastq_ftp" \
    | tail -n +2 | while IFS=$'\t' read -r run fastqs; do
        IFS=';' read -ra files <<< "$fastqs"
        for url in "${files[@]}"; do
          fname="$outdir/$(basename "$url")"
          if [[ -s "$fname" ]]; then
            echo "  - $fname exists, skipping."
          else
            echo "  - downloading $url"
            wget -c -O "$fname" "ftp://$url"
          fi
        done
      done
done
