#!/bin/bash

source $WZSEQ_ENTRY
wzref_hg38
pipeline_prepare

cd /mnt/isilon/zhou_lab/projects/20200108_biscuit_test/TruSeq_IMR90

while read sname; do
  jump_comments
  
  ## alignment
  pipeline_depend none
  fastq1=${sname}_1_L000_R1_001.fastq.gz
  fastq2=${sname}_1_L000_R2_001.fastq.gz
  output_bam=bam/${sname}.bam
  memG=10; ppn=24; queue=all.q
  pipeline_eval 11 __wgbs_biscuit_align_PE_both
  bam=bam/${sname}.bam
  memG=5; ppn=1
  pipeline_eval 12 __wzseq_index_bam

  ## markdup
  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=2; queue=all.q
  pipeline_eval 13 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=1; memG=5; ppn=1; queue=all.q
  pipeline_eval 14 __wzseq_index_bam

  ## pileup
  input_bam=bam/${sname}_markdup.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=12; memG=80; ppn=20; queue=all.q
  pipeline_depend 14
  pipeline_eval 15 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2; queue=all.q
  pipeline_eval 16 __wgbs_biscuit_QC

done <<EOM
IMR90
EOM
