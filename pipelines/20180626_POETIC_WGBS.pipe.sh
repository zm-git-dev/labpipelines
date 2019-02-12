source $WZSEQ_ENTRY
wzref_hg19_noContig
pipeline_prepare

while read sname; do
  jump_comments

  ### QC fastq ####
  fastq=fastq/${sname}_R1_001.fastq.gz
  fastq_sname=${sname}_R1
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 2 __wzseq_fastqc

  fastq=fastq/${sname}_R2_001.fastq.gz
  fastq_sname=${sname}_R2
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 3 __wzseq_fastqc

  ## alignment
  pipeline_depend none
  fastq1=fastq/${sname}_R1.fastq.gz
  fastq2=fastq/${sname}_R2.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=200; ppn=28; queue=shortq
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE_POETIC
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
  output_vcf=pileup/${sname}.vcf.gz
  hour=12; memG=80; ppn=20; queue=shortq
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

  fastq=fastq/${sname}_L005_R1_001.fastq.gz
  output_bam=bam/${sname}_SE_mate1.bam
  hour=48; memG=200; ppn=28; queue=longq;
  pipeline_depend none
  pipeline_eval 50 __wgbs_biscuit_align_SE
  bam=bam/${sname}_SE_mate1.bam
  hour=1; memG=15; ppn=1; queue=longq;
  pipeline_eval 51 __wzseq_index_bam
  
  input_bam=bam/${sname}_SE_mate1.bam
  output_sname=${sname}_SE_mate1
  hour=24; memG=50; ppn=10; queue=shortq
  pipeline_depend none
  pipeline_eval 52 __wzseq_qualimap_bamqc
  
done << EOM
Normal_0813_S19
Normal_2013_S18
Normal_2714_S20
P001-0010_cfDNA_S1
P001-0012_cfDNA_S2
P001-0014_cfDNA_S3
P001-0017_cfDNA_S4
#P001-0018_cfDNA_S5
EOM

