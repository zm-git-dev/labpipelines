source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname fqname; do
  jump_comments

  ### QC fastq ####
  input_fastq1=fastq/${fqname}_R1_001.fastq.gz
  input_fastq2=fastq/${fqname}_R2_001.fastq.gz
  output_fastq1=trimmomatic/${sname}_R1_trimmomatic.fastq.gz
  output_fastq2=trimmomatic/${sname}_R2_trimmomatic.fastq.gz
  hour=48; memG=10; ppn=10; queue=shortq
  pipeline_depend none
  pipeline_eval 1 __wzseq_trimmomatic_PE

  fastq=fastq/${fqname}_R1_001.fastq.gz
  fastq_sname=${sname^^}_R1
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 2 __wzseq_fastqc

  fastq=fastq/${fqname}_R2_001.fastq.gz
  fastq_sname=${sname^^}_R2
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 3 __wzseq_fastqc

  ## alignment
  pipeline_depend none
  fastq1=fastq/${fqname}_R1_001.fastq.gz
  fastq2=fastq/${fqname}_R2_001.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=100; ppn=10; queue=shortq
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE
  bam=bam/${sname}.bam
  hour=1; memG=5; ppn=1
  pipeline_eval 12 __wzseq_index_bam

  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=2; queue=shortq
  pipeline_eval 13 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=1; memG=5; ppn=1; queue=shortq
  pipeline_eval 14 __wzseq_index_bam

  ## pileup
  input_bam=bam/${sname}_markdup.bam
  output_vcf=bam/${sname}.new.vcf.gz
  hour=12; memG=80; ppn=10; queue=shortq
  pipeline_depend 14
  pipeline_eval 15 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=bam/${sname}.new.vcf.gz
  hour=48; memG=20; ppn=2; queue=shortq
  pipeline_eval 16 __wgbs_biscuit_QC
    
done << EOM
DMSO_d7	WGBS-DMSO-d7
EPZ_d7	WGBS-EPZ-d7
EOM

