# snDropseq_pipeline

1. Dropseq preprocessing Pipeline:
   *  dropseq_prep.sh
   *  Generate count matrix from fastq files.
2. Dropseq barcode filtering and QC Pipeline
   *  dropseq_filter_qc.py
   *  List all count matrices from DGE folder, apply mito/gene/umi filter, generate QC table for all samples.
