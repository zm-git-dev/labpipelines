#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare
cd /secondary/projects/shen/projects/2018_09_23_public_data/publication/Neurospora_crassa

while read sname srr_ids species; do

  jump_comments
  hour=48; memG=10; ppn=1
  pipeline_eval 1 __wzseq_fastq_dump_SE

  WZSEQ_BISCUIT_INDEX=~/genomes/clean/$species/biscuit/$species.fa
  fastq=fastq/${srr_ids}.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=28;
  pipeline_eval 2 __wgbs_biscuit_align_SE

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=24; memG=100; ppn=10; queue=longq
  pipeline_eval 3 __wgbs_biscuit_pileup
done << EOM
Ncrassa_WT_GSM2143335	SRR3476867	neurospora_crass_NC12
EOM


