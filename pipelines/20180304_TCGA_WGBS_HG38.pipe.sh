source $WZSEQ_ENTRY
wzref_hg38
pipeline_prepare

# cd ~/zhoulab/20200101_TCGA_WGBS
# pbsgen one 'module load python/2.7; cd ~/zhoulab/20200101_TCGA_WGBS; gdc-client download -t ~/tools/gdc-client/tokens/gdc-user-token.2020-01-01T23_57_36-05_00.txt -n 12 -m gdc_manifest.2020-01-02_TCGA2.txt' -submit -ppn 12 -name download

# awk 'NR>1{print $1}' gdc_manifest.2020-01-02_TCGA.txt | while read a; do echo $a; pbsgen one "module load python/2.7; cd ~/zhoulab/20200101_TCGA_WGBS; gdc-client download -t ~/tools/gdc-client/tokens/gdc-user-token.2020-01-01T23_57_36-05_00.txt $a" -submit -ppn 1 -name download_$a; done

# module load python/2.7
# cd ~/zhoulab/20200101_TCGA_WGBS; 
# awk 'NR>1{print $1}' gdc_manifest.2020-01-02_TCGA.txt | parallel -j 20 'a={}; echo $a; gdc-client download -t ~/tools/gdc-client/tokens/gdc-user-token.2020-01-01T23_57_36-05_00.txt $a;'
# 

# gdc-client download -t ~/tools/gdc-client/tokens/gdc-user-token.2020-01-01T23_57_36-05_00.txt 25e1cdce-1c85-4111-b2c5-696394e2d3aa

while read sname sourcebams; do

  jump_comments

  hour=48; memG=20; ppn=2; queue=all.q
  sourcebams='aggregated/'$sourcebams
  pipeline_depend none
  pipeline_eval 1 __wzseq_bam2fastq

  ## alignment
  pipeline_depend none
  fastq1=fastq/${sname}.pe1.fq.gz
  fastq2=fastq/${sname}.pe2.fq.gz
  output_bam=bam/${sname}.bam
  ppn=24; queue=all.q
  pipeline_depend none
  pipeline_eval 11 __wgbs_biscuit_align_PE
  bam=bam/${sname}.bam
  ppn=1; queue=all.q
  pipeline_eval 12 __wzseq_index_bam

  input_bam=bam/${sname}.bam
  output_bam=bam/${sname}_markdup.bam
  hour=36; memG=10; ppn=2; queue=all.q
  pipeline_eval 13 __wgbs_biscuit_markdup
  bam=bam/${sname}_markdup.bam
  hour=10; memG=5; ppn=1; queue=all.q
  pipeline_eval 14 __wzseq_index_bam

  ## pileup
  # input_bam=bam/${sname}_markdup.bam
  # skip mark duplicate
  input_bam=bam/${sname}.bam
  output_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=50; ppn=2; queue=all.q
  pipeline_depend 14
  pipeline_eval 15 __wgbs_biscuit_pileup

  input_bam=bam/${sname}_markdup.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2; queue=all.q
  pipeline_eval 16 __wgbs_biscuit_QC

  input_bam=bam/${sname}_markdup.bam
  output_sname=$sname
  hour=60; memG=200; ppn=28; queue=all.q
  pipeline_depend 14
  pipeline_eval 23 __wzseq_qualimap_bamqc
  
done << EOM
TCGA_BRCA_A04X	ResultCount_MERGING_1_NIC1254A15.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A15H	ResultCount_MERGING_1_NIC1254A18.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A07I	ResultCount_MERGING_1_NIC1254A16.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A0YG	ResultCount_MERGING_1_NIC1254A17.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_STAD_6519	ResultCount_MERGING_1_NIC1254A77.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_NA0CE	TCGA-A7-A0CE-11A-21D-A148-05_805f2c82-49d5-4994-a2a5-98dd1b832f0b_Breast_normal.realign.mdups.recal.bam
TCGA_LUSC_N2722	ResultCount_MERGING_1_NIC1254A12.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_N2689	ResultCount_MERGING_1_NIC1254A50.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_N6148	ResultCount_MERGING_1_NIC1254A70.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_STAD_N6452	ResultCount_MERGING_1_NIC1254A75.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_NA20V	ResultCount_MERGING_1_NIC1254A95.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BRCA_A0CE	TCGA-A7-A0CE-01A-11D-A148-05_4ff7626b-4c2f-40fd-9df1-10c1a488402a_Breast_tumor.realign.mdups.recal.bam
TCGA_LUSC_2722	TCGA-60-2722-01A-01D-1871-05_720f6804-059c-4516-a77d-74e8a01c739b_Lung_tumor.realign.mdups.recal.bam
TCGA_LUSC_2600	ResultCount_MERGING_1_NIC1254A13.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUSC_2695	ResultCount_MERGING_1_NIC1254A14.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_COAD_A00R	ResultCount_MERGING_1_NIC1254A45.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_0128	ResultCount_MERGING_1_NIC1254A46.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_1454	ResultCount_MERGING_1_NIC1254A47.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_3477	ResultCount_MERGING_1_NIC1254A48.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_2689	ResultCount_MERGING_1_NIC1254A49.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_READ_3593	ResultCount_MERGING_1_NIC1254A51.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_1401	ResultCount_MERGING_1_NIC1254A52.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_1460	ResultCount_MERGING_1_NIC1254A53.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_GBM_1788	ResultCount_MERGING_1_NIC1254A54.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUSC_1078	ResultCount_MERGING_1_NIC1254A67.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_4630	ResultCount_MERGING_1_NIC1254A68.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_6148	ResultCount_MERGING_1_NIC1254A69.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_6215	ResultCount_MERGING_1_NIC1254A71.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_7156	ResultCount_MERGING_1_NIC1254A72.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_LUAD_6840	ResultCount_MERGING_1_NIC1254A73.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_STAD_5730	ResultCount_MERGING_1_NIC1254A76.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_STAD_6177	ResultCount_MERGING_1_NIC1254A78.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A13J	ResultCount_MERGING_1_NIC1254A89.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A20V	ResultCount_MERGING_1_NIC1254A91.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A2HQ	ResultCount_MERGING_1_NIC1254A92.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A2LA	ResultCount_MERGING_1_NIC1254A93.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A1AG	ResultCount_MERGING_1_NIC1254A94.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_BLCA_A1AA	ResultCount_MERGING_1_NIC1254A96.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_UCEC_A05J	ResultCount_MERGING_1_NIC1254A112.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_UCEC_A0K6	ResultCount_MERGING_1_NIC1254A109.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_STAD_6452	ResultCount_MERGING_1_NIC1254A110.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_UCEC_A1CK	ResultCount_MERGING_1_NIC1254A111.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_UCEC_A0G2	ResultCount_MERGING_1_NIC1254A113.hg19_rCRSchrm.fa.realign.mdups.recal.bam
TCGA_COAD_N3518	TCGA-AA-3518-11A-01D-1518-05_d85e9d86-091c-4629-876e-d1eba7b78e83_Colon_normal.bam
TCGA_UCEC_NA1CI	TCGA-AX-A1CI-11A-11D-A17H-05_ac5058ed-1d03-4555-8f5b-a1a88eae9cc4_Endo_normal.bam
TCGA_COAD_3518	TCGA-AA-3518-01A-02D-1518-05_b9591620-97e7-4355-b18c-36fa74851106_Colon_tumor.bam
TCGA_UCEC_A1CI	TCGA-AX-A1CI-01A-11D-A17H-05_e31eb5a4-d191-4bc2-bb4d-cbd3f23bc7db_Endo_tumor.bam
EOM

