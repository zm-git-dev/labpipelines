source $SLURM_ENTRY
# wzref_mm10LambdaT4
wzref_mm10
pipeline_prepare

base=/scr1/users/zhouw3/projects/20211119_5hmC_Project/
cd $base

while read sname design srr_ids; do
  jump_comments
  pipeline_depend none
  
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  days=1; memG=100; ppn=24      # 1-2h
  pipeline_eval 1 __zlab_fasterqDumpSE_20200706

  fastq=fastq/$sname.fastq.gz
  days=2; memG=100; ppn=24      # 1-2h
  pipeline_eval 2 __zlab_trimGaloreSE_20211125

  fastq=trim/$sname/${sname}_trimmed.fq.gz
  bam_dir=bam_bismarkbt2/$sname
  days=3; memG=500; ppn=24      # 12-24h
  pipeline_eval 3 __zlab_BismarkBt2SE_20211125

  in_bam=bam_bismarkbt2/$sname/$(basename $fastq)_bismark_bt2.bam
  out_bam=bam_bismarkbt2/${sname}.bam
  days=1; memG=20; ppn=3        # 2-3h
  pipeline_eval 4 __zlab_sortIndexBam_20211126
  # Hao used the following options which all seem to be equal to default
  # nohup ${bismark_dir}/bismark --bowtie2 --fastq --multicore 4 -D 15 -R 2 -L 20 -N 0 --score_min L,0,-0.2 --gzip --bam ${genome_folder} ${output_dir}/${expt_id}_trimmed.fq.gz &> ${expt_id}_bismark.out # multicore >1 is only comptabile with --bam; --non_bs_mm: This option is only available for SAM format.
  
  input=bam_bismarkbt2/${sname}.bam
  days=2; memG=30; ppn=2        # 3h
  pipeline_eval 5 __zlab_PicardMarkdup_20211125

  in_bam=bam_PicardMdup/${sname}.bam
  days=3; memG=100; ppn=10      # 4h
  pipeline_eval 6 __zlab_BismarkMethExtract_20211125
  
done << EOM
# GSM1180315	SINGLE	SRR925933,SRR925932,SRR925938,SRR925937,SRR925936,SRR925935,SRR925934
# GSM1180316	SINGLE	SRR925940,SRR925939,SRR925945,SRR925944,SRR925943,SRR925942,SRR925941
# GSM1180317	SINGLE	SRR925947,SRR925946,SRR925953,SRR925952,SRR925951,SRR925950,SRR925949,SRR925948
# GSM1180306	SINGLE	SRR925885,SRR925884,SRR925883,SRR925882,SRR925881,SRR925880
# GSM1180307	SINGLE	SRR925891,SRR925890,SRR925889,SRR925888,SRR925887,SRR925886
# GSM1180308	SINGLE	SRR925892,SRR925893,SRR925894,SRR925895,SRR925899,SRR925897,SRR925896
# GSM1541958	SINGLE	SRR1647862,SRR1647863,SRR1647864
# GSM1541959	SINGLE	SRR1647867,SRR1647866,SRR1647865
# GSM3207185	SINGLE	SRR7368849,SRR7368848
# GSM3207186	SINGLE	SRR7368851,SRR7368850
# GSM3207181	SINGLE	SRR7368841,SRR7368842
# GSM3207183  SINGLE	SRR7368845
# GSM3207182	SINGLE	SRR7368843,SRR7368844
# GSM3207184	SINGLE	SRR7368846,SRR7368847
GSM1173794	SINGLE	SRR921893,SRR921892,SRR921891,SRR921890,SRR921889,SRR921903,SRR921902,SRR921901,SRR921900,SRR921899,SRR921898,SRR921897,SRR921896,SRR921895,SRR921894
GSM1173795	SINGLE	SRR921904,SRR921905,SRR921906,SRR921907,SRR921908,SRR921909,SRR921910,SRR921911,SRR921912,SRR921913,SRR921914,SRR921915,SRR921916,SRR921917,SRR921918
EOM
