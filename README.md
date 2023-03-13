This github is used for share the analysis code for PROMPTs.   
Chromatin RNA-seq analysis code:
```
#for paired end
bash /code/chrRNAseq_pe_star.sh hg38
#for single end
bash /code/chrRNAseq_se_star.sh hg38
```
ChIP-seq analysis code:
```
bash /code/ChIPseq_se.sh species fragment_length effective_genome_size
```
Identification of active promoters and enhancers code:
```
# for identify actively promoter
bash /code/identify_promoters.sh PolII_peak.narrowPeak all_gene_bed_file
# for identify enhancer
bash /code/identify_enhancers.sh PolII_peak.narrowPeak all_gene_bed_file
```
PROMPT-Finder:
```
bash /code/PROMPT-finder/antisense_main.sh active_gene.bed treated_mark species path_to_code
```
4sU-seq and TT-seq analysis use the same analysis pipeline:
```
bash /code/TT-seq.sh species spike_in_species
```
Transcript balance analysis code:
```
# for ChrRNA-seq
python /code/RNA_balance.py PROMPT.bed active_gene.bed scalingfactors.txt CTRL.bam treated.bam
# for TT-seq
python /code/RNA_balance_no1kb.py PROMPT.bed active_gene.bed scalingfactors.txt CTRL.bam treated.bam
# for ChIP-seq
python /code/RNAPII_balance.py PROMPT.bed active_gene.bed  CTRL.bam treated.bam 
```
U1 prediction code:
```
bash /code/predict_u1.sh genome.fa gene.bed
```
Classify gene by predicted U1 site in upatream antisense code:
```

bash /code/sort_1st_distance.py
```
