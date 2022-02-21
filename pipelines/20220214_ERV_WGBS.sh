source $SLURM_ENTRY
wzref_hg38
pipeline_prepare

base=/scr1/users/zhouw3/projects/20220214_ERV_WGBS/
cd $base

while read sname; do
  jump_comments
  pipeline_depend none
  
  fastq1=fastq/${sname}_R1.fastq.gz
  fastq2=fastq/${sname}_R2.fastq.gz
  trim_galore_dir=trim
  days=2; memG=100; ppn=6      # 1-2h
  pipeline_eval 1 __zlab_trimGalorePE_20211125

  fastq1=trim/${sname}_R1_val_1.fq.gz
  fastq2=trim/${sname}_R2_val_2.fq.gz
  outdir=bam_biscuit
  days=10; memG=110; ppn=12
  pipeline_eval 11 __zlab_biscuitAlignStrandedPE_20220216

  input=bam_biscuit/${sname}.bam
  outdir=bam_biscuit_picard
  days=2; memG=30; ppn=2        # 3h
  pipeline_eval 12 __zlab_PicardMarkdup_20220219

  in_bam=bam_biscuit_picard/${sname}.bam
  out_bam=bam_biscuit_picard/${sname}.sorted.bam
  days=1; memG=20; ppn=3        # 2-3h
  pipeline_eval 13 __zlab_sortIndexBam_20211126

  input=bam_biscuit_picard/${sname}.sorted.bam
  outdir=pileup
  days=1; memG=100; ppn=10
  pipeline_eval 14 __zlab_biscuitPileup_20220219

  input=pileup/${sname}_cg.bed.gz
  outdir=pileup
  days=1; memG=100; ppn=10
  pipeline_eval 15 __zlab_convertBigwigCG_20220219
 
  fastq1=trim/${sname}_R1_val_1.fq.gz
  fastq2=trim/${sname}_R2_val_2.fq.gz
  bam_dir=bam_bismarkbt2/${sname}
  days=21; memG=500; ppn=24
  pipeline_eval 31 __zlab_BismarkBt2StrandedPE_20220214

  input=bam_bismarkbt2/${sname}.bam
  outdir=bam_bismarkbt2_picard
  days=2; memG=30; ppn=2        # 3h
  pipeline_eval 32 __zlab_PicardMarkdup_20220219

  in_bam=bam_bismarkbt2_picard/${sname}.bam
  out_bam=bam_bismarkbt2_picard/${sname}.sorted.bam
  days=1; memG=20; ppn=3        # 2-3h
  pipeline_eval 33 __zlab_sortIndexBam_20211126
  
  in_bam=bam_bismarkbt2_picard/${sname}.sorted.bam
  outdir=bam_bismarkMeth
  days=3; memG=100; ppn=10      # 4h
  pipeline_eval 34 __zlab_BismarkMethExtract_20211125
  
  # # fastq=trim/$sname/${sname}_trimmed.fq.gz
  # # bam_dir=bam_bismarkbt2/$sname
  # # days=3; memG=500; ppn=24      # 12-24h
  # # pipeline_eval 32 __zlab_BismarkBt2SE_20211125
  
done << EOM
MCF7KD
EOM
