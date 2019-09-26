#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare

cd /secondary/projects/zhang/projects/2017_04_13_methylation_evolution/Qu_2018_GenomeRes
while read sname srr_ids species; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1
  pipeline_eval 1 __wzseq_fastq_dump_PE

  WZSEQ_BISCUIT_INDEX=/secondary/projects/jones/projects/2018_09_13_genome_assemblies/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20; queue=longq
  fastq1=fastq/${sname}_R1.fastq.gz
  fastq2=fastq/${sname}_R2.fastq.gz
  pipeline_eval 2 __wgbs_biscuit_align_PE_both
  bam=bam/${sname}.bam
  hour=5; memG=5; ppn=1
  pipeline_eval 3 __wzseq_index_bam

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=48; memG=80; ppn=10; queue=longq
  pipeline_eval 5 __wgbs_biscuit_pileup

done << EOM
Gorilla_gorilla_Sperm  SRR3289664,SRR3289665,SRR3289666,SRR3289667	gorilla_gorilla_gorGor4
Rattus_norvegicus_Sperm  SRR3289668,SRR3289669,SRR3289670,SRR3289671	rattus_norvegicus_Rnor_6.0
Canis_lupus_familiaris_Sperm_doberman  SRR3289672,SRR3289673,SRR3289674,SRR3289675,SRR3289676,SRR3289677	canis_familiaris_CanFam3.1
Canis_lupus_familiaris_Sperm_labrador  SRR3289678,SRR3289679,SRR3289680	canis_familiaris_CanFam3.1
Canis_lupus_familiaris_Sperm_portuguese_water_dog  SRR3289681,SRR3289682,SRR3289683	canis_familiaris_CanFam3.1
EOM
