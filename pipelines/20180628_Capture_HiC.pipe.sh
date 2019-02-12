source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname fastqs; do
  jump_comments

  hour=48; memG=100; ppn=20; queue=shortq;
  pipeline_depend none
  pipeline_eval 1 __hic_hicpro
  
done <<EOF
HIC2K02	HIC2K02_L000_R1_001.fastq.gz,HIC2K02_L000_R2_001.fastq.gz
HIC2WT	HIC2WT_L000_R1_001.fastq.gz,HIC2WT_L000_R2_001.fastq.gz
EOF

