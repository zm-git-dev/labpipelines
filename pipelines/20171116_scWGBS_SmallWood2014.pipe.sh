#!/bin/bash

source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

while read sname pe_or_se srr_ids; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  memG=20; ppn=5
  pipeline_eval 1 __sra_fasterq_dump_PE_20200706
done << EOM
MII_ovulated_oocytes_GSM1370522 PE      SRR1248444
MII_ovulated_oocytes_GSM1370523 PE      SRR1248445
MII_ovulated_oocytes_GSM1370524 PE      SRR1248446
MII_ovulated_oocytes_GSM1370525 PE      SRR1248447
MII_ovulated_oocytes_GSM1370526 PE      SRR1248448
MII_ovulated_oocytes_GSM1370527 PE      SRR1248449
MII_ovulated_oocytes_GSM1370528 PE      SRR1248450
MII_ovulated_oocytes_GSM1370529 PE      SRR1248451
MII_ovulated_oocytes_GSM1370530 PE      SRR1248452
MII_ovulated_oocytes_GSM1370531 PE      SRR1248453
# MII_ovulated_oocytes_GSM1370532 PE      SRR1248454
# MII_ovulated_oocytes_GSM1370533 PE      SRR1248455
# MII_ovulated_oocytes_GSM1370534 PE      SRR1248456
# mouse_ESCs_E14_GSM1370535       PE      SRR1248457
# mouse_ESCs_E14_GSM1370536       PE      SRR1248458
# mouse_ESCs_E14_GSM1370537       PE      SRR1248459
# mouse_ESCs_E14_GSM1370538       PE      SRR1248460
# mouse_ESCs_E14_GSM1370539       PE      SRR1248461
# mouse_ESCs_E14_GSM1370540       PE      SRR1248462
# mouse_ESCs_E14_GSM1370541       PE      SRR1248463
# mouse_ESCs_E14_GSM1370542       PE      SRR1248464
# mouse_ESCs_E14_GSM1370543       PE      SRR1248465
# mouse_ESCs_E14_GSM1370544       PE      SRR1248466
# mouse_ESCs_E14_GSM1370545       PE      SRR1248467
# mouse_ESCs_E14_GSM1370546       PE      SRR1248468
# mouse_ESCs_E14_GSM1370547       PE      SRR1248469
# mouse_ESCs_E14_GSM1370555       PE      SRR1248477
# mouse_ESCs_E14_GSM1370556       PE      SRR1248478
# mouse_ESCs_E14_GSM1370557       PE      SRR1248479
# mouse_ESCs_E14_GSM1370558       PE      SRR1248480
# mouse_ESCs_E14_GSM1370559       PE      SRR1248481
# mouse_ESCs_E14_GSM1370560       PE      SRR1248482
# mouse_ESCs_E14_GSM1370561       PE      SRR1248483
# mouse_ESCs_E14_GSM1370562       PE      SRR1248484
# mouse_ESCs_E14_GSM1370563       PE      SRR1248485
# mouse_ESCs_E14_GSM1370564       PE      SRR1248486
# mouse_ESCs_E14_GSM1370565       PE      SRR1248487
# mouse_ESCs_E14_GSM1370566       PE      SRR1248488
# mouse_ESCs_E14_GSM1370567       PE      SRR1248489
# mouse_ESCs_E14_GSM1370568       PE      SRR1248490
# mouse_ESCs_E14_GSM1370569       PE      SRR1248491
# mouse_ESCs_E14_GSM1370570       PE      SRR1248492
# mouse_ESCs_E14_GSM1370571       PE      SRR1248493
# FACS_Beads_GSM1370548   PE      SRR1248470
# FACS_Beads_GSM1370549   PE      SRR1248471
# Empty_ESCs_GSM1370550   PE      SRR1248472
# Empty_ESCs_GSM1370551   PE      SRR1248473
# Empty_ESCs_GSM1370552   PE      SRR1248474
# Empty_MIIs_GSM1370553   PE      SRR1248475
# Empty_MIIs_GSM1370554   PE      SRR1248476
# MII_ovulated_oocytes_deep_GSM1413827	PE	SRR1411188
# MII_ovulated_oocytes_deep_GSM1413828	PE	SRR1411189
EOM

while read sname; do
  jump_comments

  pipeline_dependlevel 1

  fastq1=fastq/${sname}_R1.fastq.gz
  fastq2=fastq/${sname}_R2.fastq.gz
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20;
  pipeline_eval 2 __wgbs_biscuit_align_PE_both

  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=1
  pipeline_eval 3 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=1; memG=5; ppn=1
  pipeline_eval 4 __wzseq_index_bam

  input_bam=bam/${sname}_markdup.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=80; ppn=10
  pipeline_eval 5 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2
  pipeline_eval 6 __wgbs_biscuit_QC

  input_bam=bam/${sname}.bam
  hour=24; memG=10; ppn=1
  pipeline_eval 8 __wzseq_uniformity_cpg
  
  input_bam=bam/${sname}.bam
  hour=24; memG=10; ppn=1
  pipeline_eval 9 __biscuit_retention

  fastq=fastq/${sname}_R1.fastq.gz
  trim_galore_dir=fastq/${sname}_trim_galore
  hour=20; memG=5; ppn=1
  pipeline_eval 12 __wzseq_trim_galore_SE

  customized_command="
  cat fastq/${sname}_trim_galore/${sname}_R1_trimmed.fq.gz fastq/${sname}_R2.fastq.gz >fastq/${sname}_trim_galore/${sname}_merged.fq.gz"
  hour=20; memG=5; ppn=1
  pipeline_eval 13 __wzseq_customize

  fastq=fastq/${sname}_trim_galore/${sname}_merged.fq.gz
  direction="--non_directional"
  bismark_bt2_dir=bam/${sname}_bismark_bt2
  bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_merged.fq.gz_bismark_bt2.bam
  bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
  hour=200; memG=180; ppn=28
  # hour=24; memG=10; ppn=1
  pipeline_eval 14 __wgbs_bismark_bowtie2_SE

  input_bam=bam/${sname}_bismark_bt2.bam
  hour=48; memG=10; ppn=1
  pipeline_eval 15 __wgbs_bismark_deduplicate

  # input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
  # hour=24; memG=10; ppn=1
  # pipeline_eval 16 __wzseq_uniformity_cpg

  input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
  hour=48; memG=10; ppn=1
  pipeline_eval 17 __wgbs_bismark_methylextraction
done <<EOM
# MII_ovulated_oocytes_GSM1370522
# MII_ovulated_oocytes_GSM1370523
# MII_ovulated_oocytes_GSM1370524
# MII_ovulated_oocytes_GSM1370525
# MII_ovulated_oocytes_GSM1370526
# MII_ovulated_oocytes_GSM1370527
# MII_ovulated_oocytes_GSM1370528
# MII_ovulated_oocytes_GSM1370529
# MII_ovulated_oocytes_GSM1370530
# MII_ovulated_oocytes_GSM1370531
# MII_ovulated_oocytes_GSM1370532
# MII_ovulated_oocytes_GSM1370533
# MII_ovulated_oocytes_deep_GSM1413827
# MII_ovulated_oocytes_deep_GSM1413828
# mouse_ESCs_E14_GSM1370535
# mouse_ESCs_E14_GSM1370536
# mouse_ESCs_E14_GSM1370537
# mouse_ESCs_E14_GSM1370538
# mouse_ESCs_E14_GSM1370539
# mouse_ESCs_E14_GSM1370540
# mouse_ESCs_E14_GSM1370541
# mouse_ESCs_E14_GSM1370542
# mouse_ESCs_E14_GSM1370543
# mouse_ESCs_E14_GSM1370544
# mouse_ESCs_E14_GSM1370545
# mouse_ESCs_E14_GSM1370546
# mouse_ESCs_E14_GSM1370555
# mouse_ESCs_E14_GSM1370556
mouse_ESCs_E14_GSM1370557
mouse_ESCs_E14_GSM1370558
mouse_ESCs_E14_GSM1370559
# mouse_ESCs_E14_GSM1370560
# mouse_ESCs_E14_GSM1370561
# mouse_ESCs_E14_GSM1370562
# mouse_ESCs_E14_GSM1370563
# mouse_ESCs_E14_GSM1370564
# mouse_ESCs_E14_GSM1370565
# mouse_ESCs_E14_GSM1370566
# mouse_ESCs_E14_GSM1370567
# mouse_ESCs_E14_GSM1370568
# mouse_ESCs_E14_GSM1370569
# mouse_ESCs_E14_GSM1370570
# mouse_ESCs_E14_GSM1370571
EOM

## the following are two bulk sequencing, we leave them out since they are very big
## MII_ovulated_oocytes_GSM1370534
## mouse_ESCs_E14_GSM1370547

