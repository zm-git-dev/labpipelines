source ~/wzlib/bash/wzseq.sh
wzref_hg38
pipeline_prepare

while read sname sourcebams; do

  jump_comments

  hour=48; memG=20; ppn=2; queue="default"
  pipeline_depend none
  pipeline_eval 1 __wzseq_bam2fastq

  ## alignment
  pipeline_depend none
  fastq1=fastq/${sname}.pe1.fq.gz
  fastq2=fastq/${sname}.pe2.fq.gz
  output_bam=bam/${sname}.bam
  hour=100; memG=200; ppn=28; queue="longq"
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE
  bam=bam/${sname}.bam
  hour=10; memG=5; ppn=1; queue=longq
  pipeline_eval 12 __wzseq_index_bam

  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=2; queue=longq
  pipeline_eval 13 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=10; memG=5; ppn=1; queue=longq
  pipeline_eval 14 __wzseq_index_bam

  ## pileup
  input_bam=bam/${sname}_markdup.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=50; ppn=5; queue=shortq
  pipeline_depend 14
  pipeline_eval 15 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2; queue=longq
  pipeline_eval 16 __wgbs_biscuit_QC

  input_bam=bam/${sname}_markdup.bam
  output_sname=$sname
  hour=60; memG=200; ppn=28; queue=longq
  pipeline_depend 14
  pipeline_eval 23 __wzseq_qualimap_bamqc
  
done << EOM
TCGA_BRCA_A04X	/primary/projects/synology/merges/2013-03-06_1323_merge_NIC1254A15/results/MERGING/MERGING_1_NIC1254A15/ResultCount_MERGING_1_NIC1254A15.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A15H	/primary/projects/synology/merges/2013-02-22_1544_merge_NIC1254A18/results/MERGING/MERGING_1_NIC1254A18/ResultCount_MERGING_1_NIC1254A18.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A07I	/primary/projects/synology/merges/2013-03-06_1325_merge_NIC1254A16/results/MERGING/MERGING_1_NIC1254A16/ResultCount_MERGING_1_NIC1254A16.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A0YG	/primary/projects/synology/merges/2013-03-06_1326_merge_NIC1254A17/results/MERGING/MERGING_1_NIC1254A17/ResultCount_MERGING_1_NIC1254A17.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_STAD_6519	/primary/projects/synology/merges/2014-02-12_1509_merge_NIC1254A77/results/MERGING/MERGING_1_NIC1254A77/ResultCount_MERGING_1_NIC1254A77.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_UCEC_A0G2	/primary/projects/synology/merges/2013-06-26_1720_merge_NIC1254A113/results/MERGING/MERGING_1_NIC1254A113/ResultCount_MERGING_1_NIC1254A113.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_COAD_N3518	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-AA-3518-11A-01D-1518-05_d85e9d86-091c-4629-876e-d1eba7b78e83_Colon_normal.realign.mdups.recal.bam
# TCGA_BRCA_NA0CE	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-A7-A0CE-11A-21D-A148-05_805f2c82-49d5-4994-a2a5-98dd1b832f0b_Breast_normal.realign.mdups.recal.bam
TCGA_UCEC_NA1CI	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-AX-A1CI-11A-11D-A17H-05_ac5058ed-1d03-4555-8f5b-a1a88eae9cc4_Endo_normal.realign.mdups.recal.bam
TCGA_LUSC_N2722	/primary/projects/synology/merges/2014-02-14_1204_merge_NIC1254A12/results/MERGING/MERGING_1_NIC1254A12/ResultCount_MERGING_1_NIC1254A12.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_N2689	/primary/projects/synology/merges/2013-06-26_1652_merge_NIC1254A50/results/MERGING/MERGING_1_NIC1254A50/ResultCount_MERGING_1_NIC1254A50.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_N6148	/primary/projects/synology/merges/2013-06-26_1752_merge_NIC1254A70/results/MERGING/MERGING_1_NIC1254A70/ResultCount_MERGING_1_NIC1254A70.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_STAD_N6452	/primary/projects/synology/merges/2014-02-13_1212_merge_NIC1254A75/results/MERGING/MERGING_1_NIC1254A75/ResultCount_MERGING_1_NIC1254A75.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_NA20V	/primary/projects/synology/merges/2013-06-26_1700_merge_NIC1254A95/results/MERGING/MERGING_1_NIC1254A95/ResultCount_MERGING_1_NIC1254A95.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_COAD_3518	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-AA-3518-01A-02D-1518-05_b9591620-97e7-4355-b18c-36fa74851106_Colon_tumor.realign.mdups.recal.bam
TCGA_BRCA_A0CE	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-A7-A0CE-01A-11D-A148-05_4ff7626b-4c2f-40fd-9df1-10c1a488402a_Breast_tumor.realign.mdups.recal.bam
# TCGA_UCEC_A1CI	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-AX-A1CI-01A-11D-A17H-05_e31eb5a4-d191-4bc2-bb4d-cbd3f23bc7db_Endo_tumor.realign.mdups.recal.bam
TCGA_LUSC_2722	/primary/projects/synology/tcgamerge/bissnpJan2013/TCGA-60-2722-01A-01D-1871-05_720f6804-059c-4516-a77d-74e8a01c739b_Lung_tumor.realign.mdups.recal.bam
TCGA_LUSC_2600	/primary/projects/synology/merges/2013-03-03_1308_merge_NIC1254A13/results/MERGING/MERGING_1_NIC1254A13/ResultCount_MERGING_1_NIC1254A13.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUSC_2695	/primary/projects/synology/merges/2013-03-03_1309_merge_NIC1254A14/results/MERGING/MERGING_1_NIC1254A14/ResultCount_MERGING_1_NIC1254A14.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_COAD_A00R	/primary/projects/synology/merges/2014-02-12_1516_merge_NIC1254A45/results/MERGING/MERGING_1_NIC1254A45/ResultCount_MERGING_1_NIC1254A45.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_0128	/primary/projects/synology/merges/2014-02-12_1511_merge_NIC1254A46/results/MERGING/MERGING_1_NIC1254A46/ResultCount_MERGING_1_NIC1254A46.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_GBM_1454	/primary/projects/synology/merges/2013-02-28_1837_merge_NIC1254A47/results/MERGING/MERGING_1_NIC1254A47/ResultCount_MERGING_1_NIC1254A47.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_GBM_3477	/primary/projects/synology/merges/2013-02-28_1835_merge_NIC1254A48/results/MERGING/MERGING_1_NIC1254A48/ResultCount_MERGING_1_NIC1254A48.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_2689	/primary/projects/synology/merges/2013-04-22_1340_merge_NIC1254A49/results/MERGING/MERGING_1_NIC1254A49/ResultCount_MERGING_1_NIC1254A49.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_3593	/primary/projects/synology/merges/2013-02-22_1649_merge_NIC1254A51/results/MERGING/MERGING_1_NIC1254A51/ResultCount_MERGING_1_NIC1254A51.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_GBM_1401	/primary/projects/synology/merges/2014-02-12_1512_merge_NIC1254A52/results/MERGING/MERGING_1_NIC1254A52/ResultCount_MERGING_1_NIC1254A52.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_GBM_1460	/primary/projects/synology/merges/2013-02-28_1835_merge_NIC1254A53/results/MERGING/MERGING_1_NIC1254A53/ResultCount_MERGING_1_NIC1254A53.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_1788	/primary/projects/synology/merges/2013-02-28_1836_merge_NIC1254A54/results/MERGING/MERGING_1_NIC1254A54/ResultCount_MERGING_1_NIC1254A54.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUSC_1078	/primary/projects/synology/merges/2013-06-26_1735_merge_NIC1254A67/results/MERGING/MERGING_1_NIC1254A67/ResultCount_MERGING_1_NIC1254A67.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_4630	/primary/projects/synology/merges/2013-06-26_1739_merge_NIC1254A68/results/MERGING/MERGING_1_NIC1254A68/ResultCount_MERGING_1_NIC1254A68.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_6148	/primary/projects/synology/merges/2013-06-26_1750_merge_NIC1254A69/results/MERGING/MERGING_1_NIC1254A69/ResultCount_MERGING_1_NIC1254A69.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_6215	/primary/projects/synology/merges/2013-06-26_1745_merge_NIC1254A71/results/MERGING/MERGING_1_NIC1254A71/ResultCount_MERGING_1_NIC1254A71.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_7156	/primary/projects/synology/merges/2013-06-26_1747_merge_NIC1254A72/results/MERGING/MERGING_1_NIC1254A72/ResultCount_MERGING_1_NIC1254A72.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_LUAD_6840	/primary/projects/synology/merges/2013-06-26_1819_merge_NIC1254A73/results/MERGING/MERGING_1_NIC1254A73/ResultCount_MERGING_1_NIC1254A73.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_STAD_5730	/primary/projects/synology/merges/2014-02-12_1437_merge_NIC1254A76/results/MERGING/MERGING_1_NIC1254A76/ResultCount_MERGING_1_NIC1254A76.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_STAD_6177	/primary/projects/synology/merges/2013-06-26_1807_merge_NIC1254A78/results/MERGING/MERGING_1_NIC1254A78/ResultCount_MERGING_1_NIC1254A78.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A13J	/primary/projects/synology/merges/2013-06-26_1656_merge_NIC1254A89/results/MERGING/MERGING_1_NIC1254A89/ResultCount_MERGING_1_NIC1254A89.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A20V	/primary/projects/synology/merges/2013-06-26_1657_merge_NIC1254A91/results/MERGING/MERGING_1_NIC1254A91/ResultCount_MERGING_1_NIC1254A91.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A2HQ	/primary/projects/synology/merges/2013-06-26_1706_merge_NIC1254A92/results/MERGING/MERGING_1_NIC1254A92/ResultCount_MERGING_1_NIC1254A92.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A2LA	/primary/projects/synology/merges/2013-06-26_1704_merge_NIC1254A93/results/MERGING/MERGING_1_NIC1254A93/ResultCount_MERGING_1_NIC1254A93.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A1AG	/primary/projects/synology/merges/2013-06-26_1702_merge_NIC1254A94/results/MERGING/MERGING_1_NIC1254A94/ResultCount_MERGING_1_NIC1254A94.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_BLCA_A1AA	/primary/projects/synology/merges/2013-06-26_1711_merge_NIC1254A96/results/MERGING/MERGING_1_NIC1254A96/ResultCount_MERGING_1_NIC1254A96.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_UCEC_A05J	/primary/projects/synology/merges/2013-06-26_1731_merge_NIC1254A112/results/MERGING/MERGING_1_NIC1254A112/ResultCount_MERGING_1_NIC1254A112.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_UCEC_A0K6	/primary/projects/synology/merges/2013-06-26_1734_merge_NIC1254A109/results/MERGING/MERGING_1_NIC1254A109/ResultCount_MERGING_1_NIC1254A109.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_STAD_6452	/primary/projects/synology/merges/2014-02-13_1214_merge_NIC1254A110/results/MERGING/MERGING_1_NIC1254A110/ResultCount_MERGING_1_NIC1254A110.hg19_rCRSchrm.fa.realign.mdups.recal.bam
# TCGA_UCEC_A1CK	/primary/projects/synology/merges/2013-06-26_1726_merge_NIC1254A111/results/MERGING/MERGING_1_NIC1254A111/ResultCount_MERGING_1_NIC1254A111.hg19_rCRSchrm.fa.realign.mdups.recal.bam
EOM

