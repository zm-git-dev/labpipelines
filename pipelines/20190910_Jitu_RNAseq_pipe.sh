source $WZSEQ_ENTRY
wzref_hg38
pipeline_prepare

while read sname fqname; do
  jump_comments

  ### QC fastq ####
  fq1=fastq/${sname}_R1.fastq.gz
  fq2=fastq/${sname}_R2.fastq.gz
  hour=48; memG=200; ppn=34; ppn_sort=6; queue=laird
  pipeline_depend none
  pipeline_eval 1 __rnaseq_STAR

done << EOM
12
13
1
23
24
25
27
2
30
33
34
35
36
37
38
39
3
40
41
44
45
46
4A
4
EOM

source ~/wzlib/bash/wzseq.sh
allbams="bam/*.bam"
# stranded="-s 1"
pairEnd="-P"
hour=48; memG=200; ppn=28; queue=shortq
# depend="-W depend=afterok:$all_aln_jobs"
echo $depend
pipeline_eval 4 __rnaseq_featureCounts

# wzseq
# wzref_hg19
# rnaseq_featureCounts
