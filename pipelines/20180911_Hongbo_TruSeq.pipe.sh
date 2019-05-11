#!/bin/bash

source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare


while read sname; do
  jump_comments
  
  fastq1=Hongbo_Data/fastq/${sname}_R1_001_val_1.fq.gz
  fastq2=Hongbo_Data/fastq/${sname}_R2_001_val_2.fq.gz
  bismark_bt2_dir=bam/${sname}_bismark_bt2
  bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_pe.bam
  direction="--non_directional"
  pbat="--pbat"
  bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
  hour=200; memG=180; ppn=28
  pipeline_eval 1 __wgbs_bismark_bowtie2_PE

  input_bam=bam/${sname}_bismark_bt2.bam
  hour=48; memG=10; ppn=1
  pipeline_eval 2 __wgbs_bismark_deduplicate

  input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
  hour=48; memG=10; ppn=1
  pipeline_eval 3 __wgbs_bismark_methylextraction
done <<EOM
FGT_HCT116_optimized_2ug_400bp
# FGT_HCT116_optimized_2ug_800bp
# FGT_HCT116_optimized_500ng_400bp
# FGT_HCT116_optimized_500ng_800bp
# VAI_HCT116_optimized_2ug_400bp
# VAI_HCT116_optimized_2ug_800bp
# VAI_HCT116_optimized_500ng_400bp
# VAI_HCT116_optimized_500ng_800bp
# FGT_HCT116_standard_500ng_400bp
# VAI_HCT116_standard_500ng_400bp
EOM
