source ~/wzlib/bash/wzseq.sh
wzref_mm10
pipeline_prepare

while read sname; do
  jump_comments

  ### QC fastq ####
  input_fastq1=fastq/${sname}_L006_R1_001.fastq.gz
  input_fastq2=fastq/${sname}_L006_R2_001.fastq.gz
  output_fastq1=fastq/${sname}_R1_trimmomatic.fastq.gz
  output_fastq2=fastq/${sname}_R2_trimmomatic.fastq.gz
  hour=48; memG=10; ppn=5; queue=shortq
  pipeline_depend none
  pipeline_eval 1 __wzseq_trimmomatic_PE

  fastq=fastq/${sname}_L006_R1_001.fastq.gz
  fastq_sname=${sname}_R1
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 2 __wzseq_fastqc

  fastq=fastq/${sname}_L006_R2_001.fastq.gz
  fastq_sname=${sname}_R2
  hour=12; memG=5; ppn=1; queue=shortq
  pipeline_depend none
  pipeline_eval 3 __wzseq_fastqc

  ## alignment
  pipeline_depend none
  fastq1=fastq/${sname}_L006_R1_001.fastq.gz
  fastq2=fastq/${sname}_L006_R2_001.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=200; ppn=28; queue=shortq
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE_Walid_lib
  bam=bam/${sname}.bam
  hour=1; memG=5; ppn=1; queue=shortq
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
  hour=12; memG=80; ppn=10; queue=shortq
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
FT-SA08501_S6
FT-SA08502_S7
FT-SA08503_S8
FT-SA08504_S9
FT-SA08505_S10
FT-SA08506_S11
#merged7to11
EOM

