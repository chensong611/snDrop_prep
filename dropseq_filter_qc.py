#!/usr/bin/python3

import os
import glob
import argparse
import subprocess
import linecache
import pandas as pd
import numpy as np

hg19_ref_bed="gencode.v24.primary_assembly.annotation.bed"
mm9_ref_bed="gencode.vM8.primary_assembly.annotation.bed"
seq_metric_script="python2.7 /media/Home_Raid1/songchen/Software/RSeQC-2.6.4_source/scripts/DropseqRnaSeqMetrics.py"

parser = argparse.ArgumentParser(description="Descriptions: Generate Dropseq QC metrics using mitochondria ratio, gene/UMI numbers as
 filter.")
parser.add_argument("--dropseq_path", dest="dropseq_path", help="Specify dropseq pipeline folder.")
parser.add_argument("--ref_path", dest="ref_path", help="Specify genome referrence folder.")
parser.add_argument("--species", dest="species", default='h', help="Specify the species the samples were collected from, default valu
e: h.")
parser.add_argument("--mito_ratio_high", default=1, type=float, dest="mito_ratio_high",
        help="Specify the maximum allowed mitochondria ratio, default value: 1.")
parser.add_argument("--gene_number_low", default=0, type=int, dest="gene_number_low",
        help="Specify the minimum allowed gene number, default value: 0.")
parser.add_argument("--gene_number_high", default=1000000, type=int, dest="gene_number_high",
        help="Specify the maximum allowed gene number, default value: 1000000")
parser.add_argument("--umi_number_low", default=0, type=int, dest="umi_number_low",
        help="Specify the minimum allowed UMI number, default value: 0.")
parser.add_argument("--umi_number_high", default=1000000, type=int, dest="umi_number_high",
        help="Specify the maximum allowed UMI number, default value: 1000000.")
parser.add_argument("--min_cell_express", default=0, type=int, dest="min_cell_express",
        help="Specify the minimum number of cells one gene should express from, default value: 0.")

args = parser.parse_args()

if args.species == 'h':
    prefix="MT-"
    ref_bed=hg19_ref_bed
elif args.species == 'm':
    prefix="mt-"
    ref_bed=mm9_ref_bed
else:
    print("Unknow species!")


print("Working under "+ os.getcwd()+" now!")
if not os.path.exists("Reports"):
    os.makedirs("Reports")
files=glob.glob("DGE/*counts.tsv")

if len(files)==0:
    print("No sample count matrix found!")
    exit(1)

QC_table=open("Reports/QCTable.csv", "w")
QC_table.write("Sample,TotalReads,polyARatio,UniqueMap,MultiMap,UnMap,Coding,UTR,Intronic,Intergenic,Nuclei#(PF),UsefulReads,UMI,Duplication\n")

for file in files:
    sample_name=file.replace("DGE/","").replace(".counts.tsv","")
    df=pd.read_table(file, index_col=0)
    barcode_list=np.array(list(df))
    all_genes=list(df.index)

    mito_genes=[ gene for gene in all_genes if gene[0:3]==prefix ]
    all_umi_number=df.sum(axis=0)
    mito_umi_number=df.loc[mito_genes].sum(axis=0)
    mito_ratio=mito_umi_number/all_umi_number
    mito_ratio=mito_ratio.fillna(0)
    mito_filter=np.array(mito_ratio<=args.mito_ratio_high)

    nonmito_genes=[ gene for gene in all_genes if gene[0:3]!=prefix ]
    nonmito_umi_number=df.loc[nonmito_genes].sum(axis=0)
    umi_filter=np.array((nonmito_umi_number>args.umi_number_low) & (nonmito_umi_number<args.umi_number_high))
    nonmito_gene_number=df.loc[nonmito_genes].astype(bool).sum(axis=0)
    gene_filter=np.array((nonmito_gene_number>args.gene_number_low) & (nonmito_gene_number<args.gene_number_high))

    filter=np.array(mito_filter & umi_filter & gene_filter)

    if np.sum(filter)>0:
        recovered_nuclei_number=np.sum(filter)
        barcode_passfilter=barcode_list[filter]
        np.savetxt("Tmp/"+sample_name+".passfilter.barcodes.txt",barcode_passfilter,fmt='%s')

        os.system(args.dropseq_path+"DigitalExpression I=Tmp/"+sample_name+".aligned.clean.bam"
            +" O=Tmp/"+sample_name+".counts.pf.tsv SUMMARY=Reports/"+sample_name+".pf.count_summary.txt"
            +" CELL_BC_FILE=Tmp/"+sample_name+".passfilter.barcodes.txt"
            +" EDIT_DISTANCE=1 OUTPUT_READS_INSTEAD=false")

        os.system(args.dropseq_path+"DigitalExpression I=Tmp/"+sample_name+".aligned.clean.bam"
            +" O=Tmp/"+sample_name+".reads.pf.tsv SUMMARY=Reports/"+sample_name+".pf.read_summary.txt"
            +" CELL_BC_FILE=Tmp/"+sample_name+".passfilter.barcodes.txt"
            +" EDIT_DISTANCE=1 OUTPUT_READS_INSTEAD=true")

        count_matrix_passfilter=pd.read_table("Tmp/"+sample_name+".counts.pf.tsv", header=0, index_col=0)
        useful_umi=count_matrix_passfilter.sum().sum()
        read_matrix_passfilter=pd.read_table("Tmp/"+sample_name+".reads.pf.tsv", header=0, index_col=0)
        useful_reads=read_matrix_passfilter.sum().sum()
        duplication=1-float(useful_umi)/useful_reads
    else:
        recovered_nuclei_number=0
        useful_umi=0
        useful_reads=0
        duplication=0

        #read distribution
    os.system(seq_metric_script+" -i Tmp/"+sample_name+".Aligned.out.sam -r "+args.ref_path+ref_bed)

    pbs_log=sample_name+".log"
    if not os.path.isfile(pbs_log):
        pbs_log=sample_name.replace("_L001","")+".log"
        if not os.path.isfile(pbs_log):
            print(sample_name+" log file not found!")
            exit(1)

    alignment_log="Reports/"+sample_name+".Log.final.out"
    if not os.path.isfile(alignment_log):
        print(sample_name+" alignment log file not found!")
        exit(1)

    genomic_distribution="Tmp/"+sample_name+".TotalReadDistribution.txt"
    if not os.path.isfile(genomic_distribution):
        print(sample_name+" read distribution file not found!")
        exit(1)

    line=linecache.getline(pbs_log,8)
    total_reads=line.split()[-1].replace(',','')
    line=linecache.getline(pbs_log,9)
    polyA_ratio=line.split('(')[1].split(')')[0]

    line=linecache.getline(alignment_log,10)
    unique_map=line.split('\t')[1].split('\n')[0]
    line=linecache.getline(alignment_log,25)
    multi_map=line.split('\t')[1].split('\n')[0]
    multi_map=float(multi_map.strip('%'))/100
    line=linecache.getline(alignment_log,27)
    many_map=line.split('\t')[1].split('\n')[0]
    many_map=float(many_map.strip('%'))/100
    line=linecache.getline(alignment_log,30)
    short_map=line.split('\t')[1].split('\n')[0]
    short_map=float(short_map.strip('%'))/100
    line=linecache.getline(alignment_log,31)
    other_map=line.split('\t')[1].split('\n')[0]
    other_map=float(other_map.strip('%'))/100

    #specific number of whitespace
    line=linecache.getline(genomic_distribution,2)
    total_tags=line.split('      ')[3].strip()
    total_tags=float(total_tags)
    line=linecache.getline(genomic_distribution,3)
    assigned_tags=line.split('      ')[1].strip()
    assigned_tags=int(assigned_tags)
    line=linecache.getline(genomic_distribution,6)
    cds=line.split('      ')[3]
    cds=int(cds)
    line=linecache.getline(genomic_distribution,7)
    utr5=line.split('      ')[3]
    utr5=int(utr5)
    line=linecache.getline(genomic_distribution,8)
    utr3=line.split('      ')[3]
    utr3=int(utr3)
    line=linecache.getline(genomic_distribution,9)
    intron=line.split('      ')[3].strip()
    intron=int(intron)
    line=linecache.getline(genomic_distribution,12)
    tss5=line.split('    ')[4]
    tss5=int(tss5)
    line=linecache.getline(genomic_distribution,15)
    tss3=line.split('    ')[3]
    tss3=int(tss3)

    record=[sample_name,total_reads,polyA_ratio,unique_map,str(multi_map+many_map),str(short_map+other_map),
            str(cds/total_tags),str((utr5+utr3)/total_tags),str(intron/total_tags),
            str((tss5+tss3+total_tags-assigned_tags)/total_tags),str(recovered_nuclei_number),
            str(useful_reads),str(useful_umi),str(duplication)+'\n']
    QC_table.write(','.join(record))


QC_table.close()
