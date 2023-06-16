# May23_Flooding_mRNAseq
ZMT-written scripts trim, align, and quantify PE50 Flooding mRNA-seq data

INPUT: demultiplexed fastq files 

OUTPUTs: 
  catenated, gzipped fastq files 
  Pre-adapter trimmed Fastqc reports
  Adapter-trimmed, gzipped fastq files
  Post-adapter trimmed Fastqc reports
  STAR transcriptome indices for transcript-model assemblies of interest
  Primary sequence alignments (.bam)
  Per-isoform and per-gene model RSEM count tables
  MultiBamSummary (heatmaps and PCAs) for within-genotype, across-treatment replicate comparison
  Differential Expression tables (from DESeq2)
