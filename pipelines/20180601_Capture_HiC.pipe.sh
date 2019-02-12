source ~/wzlib/bash/wzseq.sh
wzref_hg19
pipeline_prepare

while read sname fastqs; do
  jump_comments

  hour=48; memG=100; ppn=20; queue=shortq;
  pipeline_depend none
  pipeline_eval 1 __hic_hicpro
  
done <<EOF
DMSO_HiC	HiC-DMSO_R1_001.fastq.gz,HiC-DMSO_R2_001.fastq.gz
EPZ_HiC	HiC-EPZ_R1_001.fastq.gz,HiC-EPZ_R2_001.fastq.gz
EOF

