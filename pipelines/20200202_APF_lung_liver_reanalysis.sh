source $WZSEQ_ENTRY
wzref_hg38
pipeline_prepare

# cd /mnt/isilon/zhou_lab/projects/20200202_APF_Lung_Liver_WGBS
while read sname srr_ids; do
  jump_comments

  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1;
  pipeline_eval 1 __wzseq_fastq_dump_PE;

done <<EOF
Liver_N1    SRR2074675
Liver_T1    SRR2074677
Liver_N2    SRR2074679
Liver_T2    SRR2074681
Liver_N3    SRR2074683
Liver_T3    SRR2074685
Liver_N4    SRR2074687
Liver_T4    SRR2074689
Lung_N1    SRR2074691
Lung_T1    SRR2074693
Lung_N2    SRR2074695
Lung_T2    SRR2074697
Lung_N3    SRR2074699
Lung_T3    SRR2074701
EOF

