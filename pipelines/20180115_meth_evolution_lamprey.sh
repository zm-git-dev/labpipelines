#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare
cd /secondary/projects/shen/projects/2018_09_23_public_data/publication/Zhang2015BioRxiv

while read sname design srr_ids; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1
  pipeline_eval 1 __wzseq_fastq_dump_PE
done << EOM
SRR2457526	PAIRED	SRR2457526
SRR2457525	PAIRED	SRR2457525
SRR2442802	PAIRED	SRR2442802
EOM

while read sname fqname species; do
  jump_comments
  pipeline_dependlevel

  WZSEQ_BISCUIT_INDEX=/secondary/projects/jones/projects/2018_09_13_genome_assemblies/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20;
  fastq1=fastq/${fqname}_1.fastq.gz
  fastq2=fastq/${fqname}_2.fastq.gz
  pipeline_eval 2 __wgbs_biscuit_align_PE_both
  bam=bam/${sname}.bam
  hour=5; memG=5; ppn=1
  pipeline_eval 3 __wzseq_index_bam

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=48; memG=80; ppn=10; queue=longq
  pipeline_eval 5 __wgbs_biscuit_pileup

done << EOM
# Sperm_SRR2457526	SRR2457526	petromyzon_marinus_Pmarinus_7.0
# Heart_SRR2457525	SRR2457525	petromyzon_marinus_Pmarinus_7.0
Muscle_SRR2442802	SRR2442802	petromyzon_marinus_Pmarinus_7.0
EOM


