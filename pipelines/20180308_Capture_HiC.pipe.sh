source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname fastqs; do
  jump_comments

  hour=48; memG=100; ppn=20; queue=shortq;
  pipeline_depend none
  pipeline_eval 1 __hic_hicpro
  
done <<EOF
#K37	k37Hox_L000_R1_001.fastq.gz,k37Hox_L000_R2_001.fastq.gz
#K55	k55Hox_L000_R1_001.fastq.gz,k55Hox_L000_R2_001.fastq.gz
EPZHiC  EPZHiC_L000_R1_001.fastq.gz,EPZHiC_L000_R2_001.fastq.gz
WTHiC   WTHiC_L000_R1_001.fastq.gz,WTHiC_L000_R2_001.fastq.gz
EOF

