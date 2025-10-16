#!/bin/bash -l

#SBATCH --job-name=BLASTp
#SBATCH --nodes=1 #number of nodes requested
#SBATCH --ntasks=16 #number of threads per node
#SBATCH --export=ALL
#              d-hh:mm:ss
#SBATCH --time=3-00:00:00 # Upper time limit for job
#SBATCH --partition ac3-compute
#SBATCH --mem=32gb # How much memory you need.
#SBATCH --output=BLASTp.%j.log   # Standard output and error log
#SBATCH --mail-type=END,FAIL # Turn on mail notification
#SBATCH --mail-user=s2673271@ed.ac.uk  # Where to send mail

# Exit script on error
set -e

# Activate the conda environment for the job
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

### STEP 1) BLAST GRC proteins against ncbi protein database  ###
echo "BLASTp for GRC proteins"
blastp \
-query /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/outputs/GRC_transcripts.fasta \
-db /mnt/loki/db/ncbi_nr/nr \
-outfmt "6 std sscinames staxids stitle" \
-evalue 0.05 \
-num_threads 16 \
-out /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/outputs/GRC_transcripts_BLAST_output.tsv

# Finish the script
exit 0