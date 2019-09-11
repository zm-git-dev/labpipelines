source ~/wzlib/bash/wzseq.sh
wzref_hg19
pipeline_prepare

while read sname fqname; do
  jump_comments

  ### QC fastq ####
  fastq=fastq/$fqname
  hour=48; memG=100; ppn=20; queue=shortq
  pipeline_depend none
  pipeline_eval 1 __rnaseq_hisat2_SE

  hour=48; memG=10; ppn=2; queue=shortq
  pipeline_eval 2 __rnaseq_splitstrand_se

  hour=48; memG=10; ppn=2; queue=shortq
  pipeline_eval 3 __rnaseq_count_rmsk_stranded
  
done << EOM
K27MD10	K27MD10_L000_R1_001.fastq.gz
K27MD15	K27MD15_L000_R1_001.fastq.gz
K27MD20	K27MD20_L000_R1_001.fastq.gz
K27MDAC5	K27MDAC5_L000_R1_001.fastq.gz
K27MPBS	K27MPBS_L000_R1_001.fastq.gz
WTDAC10	WTDAC10_L000_R1_001.fastq.gz
WTDAC15	WTDAC15_L000_R1_001.fastq.gz
WTDAC20	WTDAC20_L000_R1_001.fastq.gz
WTDAC5	WTDAC5_L000_R1_001.fastq.gz
WTPBS1	WTPBS1_L000_R1_001.fastq.gz
EOM


# wzseq
# wzref_hg19
# rnaseq_featureCounts_SE_revstranded
