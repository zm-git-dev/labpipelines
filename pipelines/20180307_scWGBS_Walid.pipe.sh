source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

while read sname; do
  jump_comments

  ### QC fastq ####
  input_fastq1=fastq/${sname}_R1_001.fastq.gz
  input_fastq2=fastq/${sname}_R2_001.fastq.gz
  output_fastq1=fastq/${sname}_R1_trimmomatic.fastq.gz
  output_fastq2=fastq/${sname}_R2_trimmomatic.fastq.gz
  hour=48; memG=10; ppn=10; queue=shortq
  pipeline_depend none
  pipeline_eval 1 __wzseq_trimmomatic_PE

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
  fastq1=fastq/${sname}_R1_001.fastq.gz
  fastq2=fastq/${sname}_R2_001.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=200; ppn=28; queue=shortq
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE_Walid_lib
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
FT-SA09888_S37_L005
FT-SA09889_S38_L005
FT-SA09890_S39_L005
FT-SA09891_S40_L005
FT-SA09892_S41_L005
FT-SA09893_S42_L005
FT-SA09894_S43_L005
FT-SA09895_S44_L005
FT-SA09896_S45_L005
FT-SA09897_S46_L005
FT-SA09898_S47_L005
FT-SA09899_S48_L005
EOM

