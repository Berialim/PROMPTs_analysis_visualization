# pipeline for antisense analyse need  bam files , database.csv 
# also require expressed_list.bed generate with pol II ChIP-seq

echo "usage: bash antisense_pipeline.sh expressed_list.bed treated_mark species path_to_code"

# parameters
if [ ! -n "$4" ];then
    echo "need parameters"
    exit 0
fi
expressed=$1
treated=$2
species=$3
path_=$4

treat_bam=`awk -F "," -v vawk="$treated" '{if ($2~vawk) print $1".bam"}' database.csv`
awk '{print $4}' $expressed > expression_list.txt
# move dictionary to work place (code)
source ~/.bashrc
bam=`awk -F"," 'NR>1 {print $1".bam"}' database.csv`


# antisense annotation 
samtools merge -@ 16 treated.bam $treat_bam
samtools index -@ 16 treated.bam
bamCoverage -b treated.bam -o treated_fw.bw -bs 10 --normalizeUsing None --filterRNAstrand forward -p 16
bamCoverage -b treated.bam -o treated_rev.bw -bs 10 --normalizeUsing None --filterRNAstrand reverse -p 16
bigWigToBedGraph treated_fw.bw treated_fw.bg
bigWigToBedGraph treated_rev.bw treated_rev.bg


python $path_/sliding_.py -f treated_fw.bg -r treated_rev.bg -sp $species >>python.report
python $path_/merge_.py  -f $expressed -i result.bed -e expression_list.txt -c /reference/$species/sizes.genome >>python.report

rm temp*

# make merged 
cat $expressed merge.bed >sas.bed
python $path_/bed_to_gff.py sas.bed AST_ST.gff
mv sas.bed AST_ST.bed

featureCounts -T 16 -p -O -s 2 -a AST_ST.gff -t transcript -g transcript_id -o featureCounts.txt $bam 2>>featurecount.report

source activate R_4
Rscript $path_/antisense.R $treated $species

source ~/.bashrc

cp Differential/list300.txt expression_list.txt Differential/quantile

python $path_/database_rep.py

for i in rep_*
do
    cd $i 
    bash $path_/heatmap.sh ../Differential/quantile/ 
    bash $path_/plotProfile.sh quantile
    cd ..
done