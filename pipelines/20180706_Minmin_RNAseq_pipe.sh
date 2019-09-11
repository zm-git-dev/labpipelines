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
EVD15	EVD15_L000_R1_001.fastq.gz
EVD39	EVD39_L000_R1_001.fastq.gz
EVD55	EVD55_L000_R1_001.fastq.gz
EVD5	EVD5_L000_R1_001.fastq.gz
EVPBS	EVPBS_L000_R1_001.fastq.gz
SKDD15	SKDD15_L000_R1_001.fastq.gz
SKDD39	SKDD39_L000_R1_001.fastq.gz
SKDD55	SKDD55_L000_R1_001.fastq.gz
SKDD5	SKDD5_L000_R1_001.fastq.gz
SKDPBS	SKDPBS_L000_R1_001.fastq.gz
EOM

# wzseq
# wzref_hg19
# rnaseq_featureCounts
