This github is used for share the analysis code for PROMPTs.
System : Ubuntu(22.04.2)

Software required:
```
Trim Galore
STAR aligner v2.7.9a
SAMtools v1.7
deepTools
featureCounts v2.0.1
DESeq2
Picard
Bedtools v2.30.0
bigWigToBedGraph
```

Python3 packages for PROMPT-Finder:
```
numpy
statsmodels
pybedtools
matplotlib
argparse
```

R packages for further differential analysis and plot:
```
DESeq2
tidyverse
ggrepel
pheatmap
rstatix
ggpubr 
eulerr
```

PROMPT-Finder:
```
bash /code/PROMPT-Finder/antisense_main.sh active_gene.bed treated_mark species path_to_code
```

We also upload code for other analysis in our paper:
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
python /code/actively_transcript.py PolII_peak.narrowPeak all_gene_bed_file
# for identify enhancer
bash /code/identify_enhancers.sh PolII_peak.narrowPeak all_gene_bed_file BAMfile
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
python /code/RNAPII_balance.py PROMPT.bed active_gene.bed outfile_prefix CTRL.bam treated.bam 
```
U1 prediction code:
```
bash /code/predict_u1.sh genome.fa gene.bed
```
Classify gene by predicted U1 site in upatream antisense code:
```
# classify gene by 1st U1 site in upstream antisense
bash /code/sort_1st_distance.py
# classify gene by U1 sites number in upstream antisense 1kb of TSS
bash /code/sort_u1_number.py
```
Visualization of data through Heatmap, average line plot, boxplots, and violin plots code:
```
# plot Heatmap
bash /code/bw_heatmap.sh bedfile oudir distance bw1 bw2 ...
# plot Profile
bash /code/bw_profile.sh bedfile oudir distance bw1 bw2 ...
# plot separate strand profile
bash /code/bw_profile_sepstrand.sh bedfile outdir distance
# plot bxoplot
Rscript /code/boxplot_sas.R
Python /code/bed_boxplot.py bedfile distance_to_TSS bam1 bam2 ...
# plot violin plots same as Transcript balance analysis 
```
