source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

for file in `ls fastq/`; do
  sname=${file//.fastq.gz/};
  #sname=${file%.fastq.gz}
  pipeline_dependlevel
  #echo $sname
  
  ppn=10;
  #fastq=fastq/${sname}*.fastq.*
  #pipeline_eval 1 __fastqc20200718;
  pipeline_eval 1 __trim_galore20200806;
  
  memG=20; ppn=10
  #fastq=fastq/${sname}.fastq.gz
  fastq=TrimGalore/${sname}_trimmed.fq.gz
  pipeline_eval 2 __bwa_aln_se_20200716;

  bam="bam/"$sname".bam"
  pipeline_eval 3 __markDupAndFilter_20200719;

  species="mm"
  pipeline_eval 4 bam2bigwig;
done;

for file in `ls filtered_bam |grep ".filtered.bam$" | grep -vE 'input|MNase'`;do
  sname=${file//.filtered.bam/};
  prefix=${sname::7}
  pipeline_dependlevel
  if [ ${prefix} == "E13.5_F" ];then
    control="E13.5_F_input.filtered.bam"
  elif [ ${prefix} == "E13.5_M" ];then
    control="E13.5_M_input.filtered.bam"
  elif [ ${prefix} == "E14.5_F" ];then
    control="E14.5_F_MNase_rep1.filtered.bam"
  elif [ ${prefix} == "E14.5_M" ];then
    control="E14.5_M_MNase_rep1.filtered.bam"
  elif [ ${prefix} == "E15.5_F" ];then
    control="E15.5_F_input.filtered.bam"
  elif [ ${prefix} == "E15.5_M" ];then
    control="E15.5_M_input.filtered.bam"
  elif [ ${prefix} == "E16.5_F" ];then
    control="E16.5_F_MNase_rep1.filtered.bam"
  else
    control="E16.5_M_MNase_rep1.filtered.bam"
  fi

  pipeline_eval 5 runMacs2WithControl;

done;
