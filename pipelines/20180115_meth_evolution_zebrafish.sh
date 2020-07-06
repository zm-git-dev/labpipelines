#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare

cd /secondary/projects/shen/projects/2018_09_23_public_data/publication/Potok2013Cell
species=danio_rerio_GRCz10
while read sname se_or_pe srr_ids; do

  jump_comments
  pipeline_dependlevel

  srr_ids=${srr_ids//,/ };
  species=danio_rerio_GRCz10
  hour=48; memG=10; ppn=1
  if [[ $se_or_pe == "PAIRED" ]]; then
    pipeline_eval 1 __wzseq_fastq_dump_PE
  else
    pipeline_eval 1 __wzseq_fastq_dump_SE
  fi

  WZSEQ_BISCUIT_INDEX=~/genomes/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=28; queue=longq
  if [[ $se_or_pe == "PAIRED" ]]; then
    fastq1=fastq/${sname}_R1.fastq.gz
    fastq2=fastq/${sname}_R2.fastq.gz
    pipeline_eval 2 __wgbs_biscuit_align_PE_both
  else
    fastq=fastq/${sname}.fastq.gz
    pipeline_eval 2 __wgbs_biscuit_align_SE_both
  fi
  bam=bam/${sname}.bam
  hour=5; memG=5; ppn=1; queue=shortq
  pipeline_eval 3 __wzseq_index_bam

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=24; memG=100; ppn=10; queue=longq
  pipeline_eval 4 __wgbs_biscuit_pileup

done << EOM
SRX257168_256	SINGLE	SRR800071,SRR800072,SRR800073,SRR800074
SRX258329_256_paired	PAIRED	SRR800069,SRR800070
SRX257172_muscle	PAIRED	SRR800080,SRR800081
SRX257169_sphere	PAIRED	SRR800075,SRR800077,SRR800078
SRX257166_64	PAIRED	SRR800066,SRR800067
SRX257164_2to16	PAIRED	SRR800061,SRR800062,SRR800063,SRR800065
SRX257163_egg	PAIRED	SRR800059,SRR800060
SRX257162_sperm	PAIRED	SRR800056,SRR800057,SRR800058
EOM


