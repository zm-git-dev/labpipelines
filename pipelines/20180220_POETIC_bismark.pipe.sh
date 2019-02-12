source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname; do
  jump_comments

  ## bismark
  fastq1=fastq/${sname}_R1.fastq.gz
  fastq2=fastq/${sname}_R2.fastq.gz
  trim_galore_dir=fastq/${sname}_trim_galore
  hour=48; memG=5; ppn=1; queue=shortq
  pipeline_eval 1 __wzseq_trim_galore_PE

  fastq1=fastq/${sname}_trim_galore/${sname}_R1_val_1.fq.gz
  fastq2=fastq/${sname}_trim_galore/${sname}_R2_val_2.fq.gz
  # direction="--non_directional"
  bismark_bt2_dir=bam/${sname}_bismark_bt2
  bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_pe.bam
  bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
  hour=200; memG=180; ppn=28; queue=longq
  # hour=48; memG=20; ppn=2
  pipeline_eval 3 __wgbs_bismark_bowtie2_PE

  input_bam=bam/${sname}_bismark_bt2/${sname}_pe.bam
  hour=48; memG=50; ppn=2; queue=shortq
  library="-p"
  pipeline_eval 4 __wgbs_bismark_deduplicate

  # input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
  input_bam=bam/${sname}_bismark_bt2/${sname}_pe.deduplicated.bam
  hour=48; memG=30; ppn=2; queue=shortq
  library="--no_overlap"
  pipeline_eval 9 __wgbs_bismark_methylextraction

  input_bam=bam_biscuit/${sname}
  # __wgbs_biscuit_allele_specific_meth
  
done << EOM
# POETIC-4
# POETIC-5
# POETIC-9
# POETIC-11
# POETIC-14
# POETIC-18
# POETIC-23
# POETIC-24
# POETIC-28
# POETIC-29
# POETIC-30
# POETIC-33
# POETIC-35
# POETIC-36
# POETIC-40
# POETIC-48
# POETIC-49
# POETIC-53
# POETIC-54
# POETIC-55
# POETIC-56
# POETIC-57
# POETIC-58
# POETIC-59
# POETIC-63
# POETIC-64
# POETIC-67
# POETIC-71
# POETIC-75
# POETIC-81
# POETIC-82
POETIC-83
# POETIC-89
# POETIC-90
# POETIC-TEST
EOM

