# 01_STAR_TPM_Nmax.sh
This script maps RNA reads from t

'''
# Extract whole Embryo 0-4h: ME1, ME2, ME3, FE1, FE2, FE3
mv \
  $SCRATCH/X204SC21050493-Z01-F002/raw_data/ME{1,2,3}/*.fq.gz \
  $SCRATCH/X204SC21050493-Z01-F002/raw_data/FE{1,2,3}/*.fq.gz \
  $SCRATCH   

# Extract whole Embryo 4-8h: ML1, ML2, ML3, FL1, FL2, FL3
mv \
  $SCRATCH/X204SC21050493-Z01-F002/raw_data/ML{1,2,3}/*.fq.gz \
  $SCRATCH/X204SC21050493-Z01-F002/raw_data/FL{1,2,3}/*.fq.gz \
  $SCRATCH  

# Sync larval/pupal: Fgerm1, Fgerm2, Fgerm3, Fbody1, Fbody2, Fbody3, Mgerm1, Mgerm2, Mgerm3, Mbody1, Mbody2, Mbody3
rsync -av \
  /mnt/loki/ross/sequencing/raw/201908_Bradysia_coprophila_latelarval_earlypupa_Illumina_RNA/F*/*.fq.gz\
  /mnt/loki/ross/sequencing/raw/201908_Bradysia_coprophila_latelarval_earlypupa_Illumina_RNA/M*/*.fq.gz\
  $SCRATCH        
'''
  
  
```
echo "Trimming reads with fastp..."
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
	echo "Trimming $base"
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz && \
  	rm ${base}_1.fq.gz ${base}_2.fq.gz
done
```
