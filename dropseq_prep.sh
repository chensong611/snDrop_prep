#!/bin/sh

dropseq_dir=bin
star_dir=
picard_dir=
ref_dir=
cutadapt_dir=

species='h'
trim='T'
exonOnly='F'
cells=2000

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

read1=$sampleName'_R1_001.fastq.gz'
read2=$sampleName'_R2_001.fastq.gz'

if [ $trim = 'T' ]; then
        mv $read1 $sampleName'_R1_001.Untrimmed.fastq.gz'
        mv $read2 $sampleName'_R2_001.Untrimmed.fastq.gz'
        $cutadapt_dir/cutadapt -a TTTTTTTTT$ -m 21 -e 0.6 --discard-untrimmed -o $read1 -p $read2 $sampleName'_R1_001.Untrimmed.fastq.gz' $sampleName'_R2_001.Untrimmed.fastq.gz'
fi

echo
date

mkdir -p Reports
mkdir -p Tmp
mkdir -p DGE

#extract cell barcodes and UMI from read1, convert pair-end fastq files to bam file, and generate fastq file for alignment
java -Xmx16g -jar $picard_dir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 SAMPLE_NAME=$sampleName OUTPUT=/dev/stdout | \
java -Xmx16g -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname TMP_DIR=Tmp | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$sampleName'.cell_tag_report.txt' BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$sampleName'.molecule_tag_report.txt' BASE_RANGE=13-21 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/FilterBAM TAG_REJECT=XQ I=/dev/stdin O=/dev/stdout | \
$dropseq_dir/TrimStartingSequence INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=Reports/$sampleName'.adapter_trimming_report.txt' SEQUENCE='AAGCAGTGGTATCAACGCAGAGTGAATGGG' MISMATCHES=0 NUM_BASES=5 | \
$dropseq_dir/PolyATrimmer INPUT=/dev/stdin OUTPUT=/dev/stdout MISMATCHES=0 NUM_BASES=6 | \
tee $sampleName'.unaligned.tagged.bam' | \
java -Xmx16g -jar $picard_dir/picard.jar SamToFastq INPUT=/dev/stdin FASTQ=$sampleName'.unaligned.tagged.fastq'


#select ref genome, and perform STAR alignment
if [ $species = 'h' ]; then
        refIndex=hgSTARIndex
        refFasta=GRCh38.primary_assembly.genome.fa
        refGTF=gencode.v24.primary_assembly.annotation.gtf
else
        refIndex=mmSTARIndex
        refFasta=GRCm38.primary_assembly.genome.fa
        refGTF=gencode.vM8.primary_assembly.annotation.gtf
fi

#mapping
date
$star_dir/STAR --runThreadN 8 --genomeDir $ref_dir/$refIndex --readFilesIn $sampleName'.unaligned.tagged.fastq' --outFileNamePrefix $sampleName. --outSAMunmapped Within


#merge aligned sam file with cell barcode/UMI tagged bam file, correct barcode synthesis error, and generate digital expression matrix
if [ $exonOnly = 'F' ]; then
        TagCommand=TagReadWithGene
else
        TagCommand=TagReadWithGeneExon
fi


java -Xmx16g -jar $picard_dir/picard.jar SortSam I=$sampleName'.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=Tmp | \
java -Xmx16g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false TMP_DIR=Tmp| \
$dropseq_dir/$TagCommand I=/dev/stdin O=$sampleName'.aligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE
$dropseq_dir/DetectBeadSynthesisErrors I=$sampleName'.aligned.gene.bam' O=$sampleName'.aligned.clean.bam' OUTPUT_STATS=Reports/$sampleName'.synthesis_stats.txt' SUMMARY=Reports/$sampleName'.synthesis_stats_summary.txt' NUM_BARCODES=$cells PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
$dropseq_dir/BAMTagHistogram I=$sampleName'.aligned.clean.bam' O=DGE/$sampleName'.readsByBarcode.txt.gz' TAG=XC
$dropseq_dir/DigitalExpression I=$sampleName'.aligned.clean.bam' O=DGE/$sampleName'.counts.tsv' SUMMARY=DGE/$sampleName'.count_summary.txt' NUM_CORE_BARCODES=$cells EDIT_DISTANCE=1

if [ $trim = 'T' ]; then
        rm $read1 $read2
        mv $sampleName'_R1_001.Untrimmed.fastq.gz'  $read1
        mv $sampleName'_R2_001.Untrimmed.fastq.gz'  $read2
fi

mv $sampleName'.Log.out' Reports/
mv $sampleName'.Log.progress.out' Reports/
mv $sampleName'.SJ.out.tab' Reports/
mv $sampleName'.Log.final.out' Reports/

mv $sampleName'.aligned.clean.bam' Tmp/
mv $sampleName'.unaligned.tagged.fastq' Tmp/
mv $sampleName'.unaligned.tagged.bam' Tmp/
mv $sampleName'.Aligned.out.sam' Tmp/
mv $sampleName'.aligned.gene.bam' Tmp/
