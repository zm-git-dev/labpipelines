source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

while read sname fqname; do
  jump_comments

  ### QC fastq ####
  pipeline_depend none
  ## alignment
  fastq1=${fqname}_R1_001_val_1.fq.gz
  fastq2=${fqname}_R2_001_val_2.fq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=200; ppn=28; queue=shortq
  pipeline_eval 11 __wgbs_biscuit_align_PE_Walid_lib

  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=2; queue=shortq
  pipeline_eval 13 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=1; memG=5; ppn=1; queue=shortq
  pipeline_eval 14 __wzseq_index_bam

  ## pileup
  input_bam=bam/${sname}_markdup.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=12; memG=50; ppn=5; queue=shortq
  pipeline_depend 14
  pipeline_eval 15 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2; queue=shortq
  pipeline_eval 16 __wgbs_biscuit_QC

  input_bam=bam/${sname}_markdup.bam
  output_sname=$sname
  hour=24; memG=50; ppn=10; queue=shortq
  pipeline_depend 14
  pipeline_eval 23 __wzseq_qualimap_bamqc
  
done << EOM
# Col-3-E2-Q20 Col-3-E2-Q20/SE5250_FT-SA27061_S65_L005
# Col-3-E2-Q28 Col-3-E2-Q28/SE5250_FT-SA27061_S65_L005
# Col-4-B10-Q20 Col-4-B10-Q20/SE5249_FT-SA27018_S93_L006
# Col-4-B10-Q28 Col-4-B10-Q28/SE5249_FT-SA27018_S93_L006
# SI-2-G1-Q20 SI-2-G1-Q20/SE5229_FT-SA24002_S72_L008
Col-3-E2-Q28 Col-3-E2-Q28/SE5250_FT-SA27061_S65_L005
Col-4-B10-Q20 Col-4-B10-Q20/SE5249_FT-SA27018_S93_L006
EOM
