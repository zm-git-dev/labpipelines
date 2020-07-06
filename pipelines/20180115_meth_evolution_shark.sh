# cd /secondary/projects/zhang/projects/2017_04_13_methylation_evolution/ElephantShark_C_milii_WGBS_GSE96683
# [sra]

#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare
cd /secondary/projects/zhang/projects/2017_04_13_methylation_evolution/ElephantShark_C_milii_WGBS_GSE96683

while read sname srr_ids; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1
  pipeline_eval 1 __wzseq_fastq_dump_PE
done << EOM
Shark_Liver_female  SRR5349298
Shark_Liver_male  SRR5349299
EOM

while read sname species; do
  jump_comments
  pipeline_dependlevel

  WZSEQ_BISCUIT_INDEX=/secondary/projects/jones/projects/2018_09_13_genome_assemblies/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20;
  fastq=fastq/${sname}.fastq.gz
  pipeline_eval 2 __wgbs_biscuit_align_SE_both

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=48; memG=80; ppn=10; queue=longq;
  pipeline_eval 5 __wgbs_biscuit_pileup

done << EOM
Shark_Liver_female	callorhinchus_milii_calMil1
Shark_Liver_male	callorhinchus_milii_calMil1
EOM


