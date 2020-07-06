#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare

cd /secondary/projects/jones/projects/2018_09_13_genome_assemblies/raw/Collaboration/european_amphioxus
while read sname design srr_ids; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1
  pipeline_eval 11 __wzseq_fastq_dump_SE
done << EOM
GSM2728827_8hpf	SINGLE	SRR5889231,SRR5889232,SRR5889233,SRR5889234,SRR5889235,SRR5889236,SRR5889237,SRR5889238
GSM2728828_15hpf	SINGLE	SRR5889239,SRR5889240,SRR5889241,SRR5889242
GSM2728829_36hpf	SINGLE	SRR5889243,SRR5889244,SRR5889245,SRR5889246,SRR5889247
GSM2728830_liver	SINGLE	SRR5889248,SRR5889249,SRR5889250,SRR5889251,SRR5889252
EOM


while read sname species; do
  jump_comments
  pipeline_dependlevel

  WZSEQ_BISCUIT_INDEX=/secondary/projects/jones/projects/2018_09_13_genome_assemblies/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20;
  fastq=fastq/${sname}.fastq.gz
  pipeline_eval 12 __wgbs_biscuit_align_SE_both
  bam=bam/${sname}.bam
  hour=5; memG=5; ppn=1
  pipeline_eval 13 __wzseq_index_bam

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=48; memG=80; ppn=10
  pipeline_eval 15 __wgbs_biscuit_pileup

done << EOM
GSM2728827_8hpf	amphioxus_Bl71nemr
GSM2728828_15hpf	amphioxus_Bl71nemr
GSM2728829_36hpf	amphioxus_Bl71nemr
GSM2728830_liver	amphioxus_Bl71nemr
EOM



