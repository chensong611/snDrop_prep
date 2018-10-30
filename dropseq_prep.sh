#!/bin/sh

dropseq_dir=/home/sochen/oasis/Software/bin
star_dir=/home/sochen/oasis/Software/STAR-2.5.2b
picard_dir=/home/sochen/oasis/Software/picard-tools-2.1.0
ref_dir=/home/sochen/oasis/Ref
cutadapt_dir=/home/sochen/oasis/Software/python/bin
kneeplot_dir=/home/sochen/oasis/Scripts

export JAVA_HOME=/home/sochen/oasis/Software/jre1.8.0_121/bin
export PATH=/home/sochen/oasis/Software/anaconda2/bin:/home/sochen/oasis/Scripts:/home/sochen/oasis/Software/bin:$JAVA_HOME:$PATH
module load python
export PYTHONPATH=/home/sochen/oasis/Software/python/lib/python2.7/site-packages:$PYTHONPATH
export LIBRARY_PATH=/home/sochen/oasis/Software/gcc-5.3.0/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/sochen/oasis/Software/gcc-5.3.0/lib64:$LD_LIBRARY_PATH
module load R

species='h'
exonOnly='F'
cells=2000
trim='T'

while getopts ":n:s:e:c:t:" options; do
    case $options in
        n ) sampleName=$OPTARG;;
        s ) species=$OPTARG;;
        e ) exonOnly=$OPTARG;;
        c ) cells=$OPTARG;;
        t ) trim=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

if [ $trim = 'T' ]; then
        mv $sampleName'_R1_001.fastq.gz' $sampleName'_R1_001.Untrimmed.fastq.gz'
        mv $sampleName'_R2_001.fastq.gz' $sampleName'_R2_001.Untrimmed.fastq.gz'
        $cutadapt_dir/cutadapt -a TTTTTTTTT$ -m 21 -e 0.6 --discard-untrimmed -o $sampleName'_R1_001.fastq.gz' -p $sampleName'_R2_001.fastq.gz' $sampleName'_R1_001.Untrimmed.fastq.gz' $sampleName'_R2_001.Untrimmed.fastq.gz'
fi

echo
date

read1=$sampleName'_R1_001.fastq.gz'
read2=$sampleName'_R2_001.fastq.gz'

mkdir -p Reports

#extract cell barcodes and UMI from read1, convert pair-end fastq files to bam file, and generate fastq file for alignment
java -Xmx16g -jar $picard_dir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 SAMPLE_NAME=$sampleName OUTPUT=/dev/stdout | \
java -Xmx16g -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname TMP_DIR=`pwd`/Tmp | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.cell_tag_report.txt' BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.molecule_tag_report.txt' BASE_RANGE=13-21 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/FilterBAM TAG_REJECT=XQ I=/dev/stdin O=/dev/stdout | \
$dropseq_dir/TrimStartingSequence INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=`pwd`/Reports/$sampleName'.adapter_trimming_report.txt' SEQUENCE='AAGCAGTGGTATCAACGCAGAGTGAATGGG' MISMATCHES=0 NUM_BASES=5 | \
$dropseq_dir/PolyATrimmer INPUT=/dev/stdin OUTPUT=/dev/stdout MISMATCHES=0 NUM_BASES=6 | \
tee $sampleName'.unaligned.tagged.bam' | \
java -Xmx16g -jar $picard_dir/picard.jar SamToFastq INPUT=/dev/stdin FASTQ=`pwd`/$sampleName'.unaligned.tagged.fastq'


#select ref genome, and perform STAR alignment
if [ $species = 'h' ]; then
        refIndex=humGenomeIndex_99
        refFasta=GRCh38.primary_assembly.genome.fa
                refGTF=gencode.v24.primary_assembly.annotation.gtf
else
        refIndex=musGenomeIndex
        refFasta=GRCm38.primary_assembly.genome.fa
                refGTF=gencode.vM8.primary_assembly.annotation.gtf
fi


date
$star_dir/STAR --runThreadN 8 --genomeDir $ref_dir/$refIndex --readFilesIn `pwd`/$sampleName'.unaligned.tagged.fastq' --outFileNamePrefix $sampleName. --outSAMunmapped Within


#merge aligned sam file with cell barcode/UMI tagged bam file, correct barcode synthesis error, and generate digital expression matrix
if [ $exonOnly = 'F' ]; then
        TagCommand=TagReadWithGene
else
        TagCommand=TagReadWithGeneExon
fi

mkdir -p Tmp
mkdir -p DGE
mkdir -p QC

java -Xmx16g -jar $picard_dir/picard.jar SortSam I=$sampleName'.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=`pwd`/Tmp | \
java -Xmx16g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false TMP_DIR=`pwd`/Tmp| \
$dropseq_dir/$TagCommand I=/dev/stdin O=$sampleName'.aligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE
$dropseq_dir/DetectBeadSynthesisErrors I=$sampleName'.aligned.gene.bam' O=$sampleName'.aligned.clean.bam' OUTPUT_STATS=`pwd`/Reports/$sampleName'.synthesis_stats.txt' SUMMARY=`pwd`/Reports/$sampleName'.synthesis_stats_summary.txt' NUM_BARCODES=$cells PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
$dropseq_dir/BAMTagHistogram I=$sampleName'.aligned.clean.bam' O=`pwd`/DGE/$sampleName'.readsByBarcode.txt.gz' TAG=XC
$dropseq_dir/DigitalExpression I=$sampleName'.aligned.clean.bam' O=`pwd`/DGE/$sampleName'.counts.tsv' SUMMARY=`pwd`/DGE/$sampleName'.count_summary.txt' NUM_CORE_BARCODES=$cells EDIT_DISTANCE=1

if [ $trim = 'T' ]; then
        rm $sampleName'_R1_001.fastq.gz' $sampleName'_R2_001.fastq.gz'
        mv $sampleName'_R1_001.Untrimmed.fastq.gz'      $sampleName'_R1_001.fastq.gz'
        mv $sampleName'_R2_001.Untrimmed.fastq.gz'      $sampleName'_R2_001.fastq.gz'
fi

mv $sampleName'.Log.out' Reports/
mv $sampleName'.Log.progress.out' Reports/
mv $sampleName'.SJ.out.tab' Reports/

mv $sampleName'.Log.final.out' QC/

mv $sampleName'.aligned.clean.bam' Tmp/
mv $sampleName'.unaligned.tagged.fastq' Tmp/
mv $sampleName'.unaligned.tagged.bam' Tmp/
mv $sampleName'.Aligned.out.sam' Tmp/
mv $sampleName'.aligned.gene.bam' Tmp/

cd DGE
$kneeplot_dir/KneePlot.R $sampleName $(($cells*5))
