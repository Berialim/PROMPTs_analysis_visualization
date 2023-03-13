# RNAseq single end
echo "usage:"
echo "bash chrrnaseq_se_star.sh species(hg38. ce11)"
if [ ! -n "$1" ];then
	echo "no enough parameters"
	exit 8
fi

if [ ! -f "database.csv" ]; then
    printf "need database.csv for sample information, \nform:\n"
    printf "name,condition\nnam1,condition1\nname2,condition2\n...,...\n"
    exit 8
fi

species=$1
feature=transcript
name=gene_id

source /home/fanlab/.bashrc
# need all fastq file in fq files
cd fq 
for i in *.fq.gz
do
    ii=`echo $i | sed "s/.*_L0._//g"`
    # remove the sequencing mark (FXXXX_L01_sample.fq.gz)
    mv $i $ii
done
# mv ambiguous reads and undecoded reads into useless folder
mkdir useless
mv ambiguous.fq.gz useless
mv undecoded.fq.gz useless
x=`ls *_1.fq.gz | sed 's/_1.fq.gz//g'`
cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore --paired -q 30 fq/${line}_1.fq.gz  fq/${line}_2.fq.gz  -o trimmed &
done
wait
# remove the mark of trimmed 
cd trimmed
for i in *.fq.gz
do
    mv ${i} ${i%%_val*}.fq.gz
done
cd ..
# remove the rRNA reads
mkdir de_rRNA
touch rRNA_report
for line in $x
do
    echo 'rRNA_sample:'$line >> rRNA_report
    bowtie -p 16 --un de_rRNA/${line}.fq -x /reference/${species}/ribosome/rDNA -1 trimmed/${line}_1.fq.gz -2 trimmed/${line}_2.fq.gz -S useless.sam 2>>rRNA_report
done
# plot rRNA percentage 
python /home/fanlab/code/rRNA_plot.py
rm *.sam
# get the quality of fastq
mkdir qc
mkdir fastqc
for line in $x
do
    fastqc de_rRNA/${line}_1.fq de_rRNA/${line}_2.fq -o fastqc &
done
wait
mkdir aligned
for line in $x
do
    STAR --genomeDir /reference/hg38/star_index/ --runThreadN 20 \
    --readFilesIn de_rRNA/${line}_1.fq de_rRNA/${line}_2.fq  --outFileNamePrefix aligned/${line} --outSAMtype SAM \
    --outSAMattributes Standard
    mv aligned/${line}Aligned.out.sam aligned/${line}.sam
    mv aligned/${line}Log* fastqc
    samtools sort -@ 16 aligned/${line}.sam -o aligned/${line}.bam
    rm aligned/${line}.sam
    samtools index -@ 16 aligned/${line}.bam
done

mkdir post_QC
mkdir tempQC
# RSeqc quality control for bam file 
for line in $x
do
     bam_stat.py -i aligned/${line}.bam >> post_QC/${line}_state.txt&
     read_distribution.py -r /reference/$species/RefSeq.bed  -i aligned/${line}.bam  >> post_QC/${line}_read_distribution.txt&
wait
done
rm -r tempQC
mv post_QC/* fastqc
multiqc fastqc 
cd aligned
cp ../*.csv ./
mkdir ../bw
mkdir ../bw_nosep
# generate bigwig files
for i in *.bam
do
    bamCoverage -b $i -o ../bw/${i%.*}_fw.bw --normalizeUsing RPKM --binSize 10 --filterRNAstrand forward -p 16
    bamCoverage -b $i -o ../bw/${i%.*}_rev.bw --normalizeUsing RPKM --binSize 10 --filterRNAstrand reverse -p 16
    # generate no separate bw file 
    bamCoverage -b $i -o ../bw_nosep/${i%.*}.bw --normalizeUsing RPKM --binSize 10 -p 16
done
cd ../bw
for i in *_rev.bw
do
    bash /home/fanlab/code/neg_bw.sh $i &
done
cd ..

cp database.csv bw
cp *.csv aligned
cd aligned

csv=`ls *.csv|sed "s/\.csv//g"`
for i in *.csv
do
    mkdir ${i%.*}
    mv $i ${i%.*}/database.csv
done
for c in $csv
do
    cd $c
    # generate featureCounts.txt for DESeq2 analysis
    bam=` awk -F"," 'NR>1{print "../"$1".bam"}' database.csv`
    featureCounts -T -p 16 -O -s 2 -a /reference/$species/ref_nmr.gtf -t $feature -g $name -o featureCounts.txt $bam
    # analyze with DESeq2 in R
    conda activate R_4
    Rscript /home/fanlab/code/RNAseq_gene_id.R
    cd ..
    source /home/fanlab/.bashrc
done
