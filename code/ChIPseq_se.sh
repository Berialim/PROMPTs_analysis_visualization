##chipseq_se
mkdir aligned
mkdir trimmed
mkdir fastqc1

# move all fastq in fq dictionary
echo "usage:"
echo "bash chipseq_se.sh species(e.g. hg38) fragment_length(from QSep) effectiveGenomeSize(for deeptools bamCoverage)"
if [ ! -n "$3" ];then
        echo "not enough parameter"
        exit 8
fi
species=$1
fragment_length=$2
genomesize=$3

source /home/fanlab/.bashrc
cd fq
# rename raw data with prefix FXXXX_L0*_X.fq.gz
for i in *.fq.gz
do
    mv $i ${i##*L0*_}
done
mkdir useless
mv ambiguous.fq.gz useless
mv undecoded.fq.gz useless
x=`ls *.fq.gz | sed 's/.fq.gz//g'`
cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore -q 30 fq/${line}.fq.gz  -o trimmed &
done
wait
cd trimmed
for i in *.fq.gz
do
    mv ${i} ${i%_*}.fq.gz
done
cd ..

for i in trimmed/*.fq.gz
do 
    fastqc ${i} -o fastqc1& 
done
wait
multiqc fastqc1 -o qc

cd trimmed
for i in *.fq.gz
do
    echo "sample: ${i%%.*}" >> ../qc/align_report
    # align to genome
    bowtie2 -p 16 -x /reference/$species/bowtie2/genome -U ${i}  -S ../aligned/${i%%.*}.sam 2>> ../qc/align_report
    # generate bam file from sam file 
    samtools sort -@ 16 -o ../aligned/${i%%.*}.bam  ../aligned/${i%%.*}.sam
    # remove duplicates
    picard MarkDuplicates I=../aligned/${i%%.*}.bam  O=../aligned/${i%%.*}.markdup.bam  M=../aligned/${i%%.*}.markdup.txt
    # remove temp files
    rm ../aligned/${i%%.*}.sam
    rm ../aligned/${i%%.*}.bam
    mv ../aligned/${i%%.*}.markdup.bam ../aligned/${i%%.*}.bam
    # generate bam index
    samtools index -@ 16 ../aligned/${i%%.*}.bam
done
cd ../qc
python /home/fanlab/code/alignment_se.py

# generatre bigwig file for read coverage in genome(visualize in igv)
mkdir ../bw
cd ../aligned
for i in *.bam
do
    bamCoverage --bam ${i} -o ../bw/${i%.*}.bw --binSize 10 --normalizeUsing RPGC \
     --effectiveGenomeSize $genomesize -e $fragment_length -p 16
done

rm -r trimmed


