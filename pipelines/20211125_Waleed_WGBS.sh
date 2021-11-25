source $SLURM_ENTRY
wzref_mm10
pipeline_prepare

base=/scr1/users/zhouw3/projects/20211119_5hmC_Project/

while read sname; do
  jump_comments

  hour=48; memG=64; ppn=24; queue=defq
  pipeline_depend none

  trim_galore_dir=fastq/trimmed/
  pipeline_eval 1 __zlab_trimGaloreSE_20211125

  fastq=trim/${sname}
  direction="--non_directional"
  bismark_bt2_dir=bam/${sname}
  bismark_bt2_bam_unsorted=bam_bismarkbt2/${sname}/${sname}_merged.fq.gz_bismark_bt2.bam
  bismark_bt2_bam_final=bam/${sname}.bam
  hour=48; memG=180; ppn=10; queue=defq
  pipeline_eval 2 __zlab_BismarkBt2SE_20211125

  #   cat <<'EOF' | sbatch --mem=64G -c 10 -t 3-2
  # #!/bin/bash
  # cd $base
  # input=/mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/TAB-seq/GSM1180307/bam/GSM1180307_bis_bt2.bam
  # java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar ~/zhoulab/labsoftware/picard/picard-2.23.2.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=picard/GSM1180307_bis_bt2.mdup.stats READ_NAME_REGEX=null INPUT=${input} OUTPUT=picard/GSM1180307_bis_bt2.bam TMP_DIR=tmp
  # samtools flagstat picard/GSM1180307_bis_bt2.bam > picard/GSM1180307_bis_bt2.bam.flagstat
  # EOF

  #   for filename in *.bam; do
  #   echo "extracting methylation for: ${filename}"
  #   /mnt/isilon/zhoulab/labsoftware/bismark/bismark_v0.14.5/bismark_methylation_extractor --no_overlap --multicore 5 --report ${filename}
  # done

done << EOM
# GSM1180315
GSM1180316
# GSM1180317
# GSM1180306
# GSM1180307
# GSM1180308
# GSM1541958
# GSM1541959
# GSM3207185
# GSM3207186
EOM


# cat ~/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_Neurons/ACE-seq/GSM3207185/SRR7368848.fastq ~/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_Neurons/ACE-seq/GSM3207185/SRR7368849.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM3207185.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_Neurons/ACE-seq/GSM3207186/SRR7368850.fastq /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_Neurons/ACE-seq/GSM3207186/SRR7368851.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM3207186.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/TAB-seq/GSM1180307/GSM1180307.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM1180307.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/TAB-seq/GSM1180308/GSM1180308.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM1180308.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/WGBS/GSM1180315/GSM1180315.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM1180315.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/WGBS/GSM1180316/GSM1180316.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM1180316.fastq.gz
# cat /mnt/isilon/zhoulab/labprojects/20210827_5hmCProject_WI/20210802_ESCs/WGBS/GSM1180317/GSM1180317.fastq | pigz -p 24 -c >~/scr1_zhouw3/projects/20211119_5hmC_Project/fastq/GSM1180317.fastq.gz
